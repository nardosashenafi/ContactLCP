import ChainRulesCore

function lcpGradForwardSysParam(lcp, x, θm, stateNum; Δt = 0.001)

    model = JuMP.Model(() -> DiffOpt.diff_optimizer(Mosek.Optimizer))
    set_silent(model)

    gn, γn, γt, M, h, Wn, Wt = ContactLCP.setSysAttributes(lcp, x)
    ContactLCP.checkContact(lcp, gn)
    ContactLCP.trimAttributes(lcp, gn, γn, γt, M, h, Wn, Wt)

    #change the mass, update the mass matrix and the gravitational force
    M = [θm 0.0; 0.0 θm]
    Minv = inv(M)
    h = [0.0, -θm*sys.g]

    E = Matrix{Float64}(I, lcp.current_contact_num, lcp.current_contact_num)

    A = lcp.Wn'*(M\lcp.Wn)
    b = lcp.Wn'*(M \ h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn

    #TODO: use env macro to directly receive state info from the system struct
    qM = x[1:stateNum]
    uA = x[stateNum+1:2stateNum]

    @variable(model, λ[1:lcp.total_contact_num] >= 0.0)
    @variable(model, v[1:stateNum])
    @variable(model, q[1:stateNum])

    @expression(model, contactForces, [lcp.sys.contactIndex[i] == 1 ? λ[i] : 0.0 for i in 1:lcp.total_contact_num] )
    @constraint(model,
        cons1,
        A*λ[lcp.sys.contactIndex .== 1] .+ lcp.Wn'*(Minv * h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn .>= 0.0
    )
    @constraint(model,
        cons2[j in 1:stateNum],
        sum(M[j,i]*(v[i] - uA[i]) - h[j]*Δt for i in 1:stateNum) .== sum(Wn[j, k]*contactForces[k] for k in 1:lcp.total_contact_num)
    )

    @constraint(model, 
        cons3, 
        q .== qM + 0.5*Δt*v
    )

    @objective(model, Min, dot(A*λ[lcp.sys.contactIndex .== 1] .+ b, λ[lcp.sys.contactIndex .== 1]))

    optimize!(model) 
    
    q_new = JuMP.value.(q)
    v_new = JuMP.value.(v)
    λ_new = JuMP.value.(λ)

    MOI.set.(  
        model, 
        DiffOpt.ForwardConstraintFunction(),    #to take the derivative wrt to θm
        cons2[2], 
        v[2] - uA[2]      #∂cons2[1]/∂θm
    )
    
    DiffOpt.forward_differentiate!(model)

    ∂λ∂θ = MOI.get.(   
            model,
            DiffOpt.ForwardVariablePrimal(),
            [model[:q]; model[:v]; model[:λ]])

    return q_new, v_new, λ_new, ∂λ∂θ
end

function lcpSolve(lcp, x, θm, stateNum;  model = JuMP.Model(Mosek.Optimizer), Δt = 0.001)

    set_silent(model)

    gn, γn, γt, M, h, Wn, Wt = ContactLCP.setSysAttributes(lcp, x)
    ContactLCP.checkContact(lcp, gn)
    ContactLCP.trimAttributes(lcp, gn, γn, γt, M, h, Wn, Wt)

    #change the mass, update the mass matrix and the gravitational force
    M = [θm 0.0; 0.0 θm]
    Minv = inv(M)
    h = [0.0, -θm*sys.g]

    E = Matrix{Float64}(I, lcp.current_contact_num, lcp.current_contact_num)

    A = lcp.Wn'*(M\lcp.Wn)
    b = lcp.Wn'*(M \ h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn

    #TODO: use env macro to directly receive state info from the system struct
    qM = x[1:stateNum]
    uA = x[stateNum+1:2stateNum]

    @variable(model, λ[1:lcp.total_contact_num] >= 0.0)
    @variable(model, v[1:stateNum])
    @variable(model, q[1:stateNum])

    @expression(model, contactForces, [lcp.sys.contactIndex[i] == 1 ? λ[i] : 0.0 for i in 1:lcp.total_contact_num] )
    @constraint(model,
        cons1,
        A*λ[lcp.sys.contactIndex .== 1] .+ lcp.Wn'*(Minv * h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn .>= 0.0
    )
    @constraint(model,
        cons2[j in 1:stateNum],
        sum(M[j,i]*(v[i] - uA[i]) - h[j]*Δt for i in 1:stateNum) .== sum(Wn[j, k]*contactForces[k] for k in 1:lcp.total_contact_num)
    )

    @constraint(model, 
        cons3, 
        q .== qM + 0.5*Δt*v
    )

    @objective(model, Min, dot(A*λ[lcp.sys.contactIndex .== 1] .+ b, λ[lcp.sys.contactIndex .== 1]))

    optimize!(model) 

   return [JuMP.value.(q), JuMP.value.(v), JuMP.value.(λ)]
end

function ChainRulesCore.frule((uA, Δx), 
                            ::typeof(lcpSolve),
                            lcp, x, θm, stateNum
                        )

    model = JuMP.Model(() -> DiffOpt.diff_optimizer(Mosek.Optimizer))
    sol = lcpSolve(lcp, x, θm, stateNum, model=model)
    q_new, v_new, λ_new = sol

    cons2 = model[:cons2]

    MOI.set.(  
        model, 
        DiffOpt.ForwardConstraintFunction(),    #to take the derivative wrt to θm
        cons2[2], 
        v_new[2] - uA[2]      #∂cons2[1]/∂θm
    )

    DiffOpt.forward_differentiate!(JuMP.backend(model))
    ∂sol∂θm =  MOI.get.(   
            model,
            DiffOpt.ForwardVariablePrimal(),
            [model[:q]; model[:v]; model[:λ]])

    return (sol, ∂sol∂θm)
end

function solveM(lcp, x0; Δt = 0.001, totalTimeStep = 1000)

    S       = Array{Array{T, 1}, 1}()
    Λ      = Array{Array{T, 1}, 1}()
    θm_actual = 0.2
    x1 = deepcopy(x0)
    #generate trajectory
    for i in 1:totalTimeStep
        expert_sol = lcpSolve(lcp, x1, θm_actual, 2)
        x1 = vcat(expert_sol[1], expert_sol[2])
        push!(S, deepcopy(x1))
        push!(Λ, deepcopy(expert_sol[3]))
    end

    θm0 = 0.9

    loss(θm) = computeLoss(lcp, S, θm, x0; totalTimeStep=totalTimeStep)
    grad_m = ForwardDiff.derivative(loss, θm0)

end

function computeLoss(lcp, expert_sol, θm, x0; totalTimeStep = 1000)

    loss = 0.0
    x1 = deepcopy(x0)
    for i in 1:totalTimeStep
        sol = lcpSolve(lcp, x1, θm, 2)
        x1 = vcat(sol[1], sol[2])
        loss += norm(expert_sol[i] - x1)
    end
    
    return loss

end