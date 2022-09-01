import ChainRulesCore, ForwardDiff
using JuMP, Mosek, MosekTools, DiffOpt, LinearAlgebra

function lcpSolve(::Type{T}, lcp, x, θm, stateNum;  model = JuMP.Model(Mosek.Optimizer), Δt = 0.001) where T

    set_silent(model)

    gn, γn, γt, M, h, Wn, Wt = ContactLCP.setSysAttributes(lcp, x)
    ContactLCP.checkContact(lcp, gn)
    ContactLCP.trimAttributes(lcp, gn, γn, γt, M, h, Wn, Wt)

    #change the mass, update the mass matrix and the gravitational force
    M = T.([θm 0.0; 0.0 θm])
    Minv = inv(M)
    h = T.([0.0, -θm*lcp.sys.g])

    E = Matrix{T}(I, lcp.current_contact_num, lcp.current_contact_num)

    A = T.(lcp.Wn'*(Minv*lcp.Wn))
    b = T.(lcp.Wn'*(Minv*h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn)

    #TODO: use env macro to directly receive state info from the system struct
    qM = T.(x[1:stateNum])
    uA = T.(x[stateNum+1:2stateNum])

    @variable(model, λ[1:lcp.total_contact_num] >= T(0.0))
    @variable(model, v[1:stateNum])
    @variable(model, q[1:stateNum])

    @expression(model, contactForces, [lcp.sys.contactIndex[i] == 1 ? λ[i] : T(0.0) for i in 1:lcp.total_contact_num] )
    @constraint(model,
        cons1,
        A*λ[lcp.sys.contactIndex .== 1] .+ b .>= T(0.0)
    )

    # @constraint(model,
    #     cons2[j in 1:stateNum],
    #     sum(M[j,i]*(v[i] - uA[i]) for i in 1:stateNum) - h[j]*Δt .== sum(Wn[j, k]*contactForces[k] for k in 1:lcp.total_contact_num)
    # )

    @constraint(model,
        cons2,
        M*(v - uA) - h*Δt .== Wn*contactForces
    )

    @constraint(model, 
        cons3, 
        q .== qM + 0.5*Δt*v
    )

    @objective(model, Min, dot(A*λ[lcp.sys.contactIndex .== 1] .+ b, λ[lcp.sys.contactIndex .== 1]))

    optimize!(model) 

   return vcat(JuMP.value.(q), JuMP.value.(v), JuMP.value.(λ)), A, uA
end

function ChainRulesCore.frule((_, _), 
                            ::typeof(lcpSolve),
                            ::Type{T}, lcp, x, θm, stateNum, Δt
                        ) where T

    model = JuMP.Model(() -> DiffOpt.diff_optimizer(Mosek.Optimizer))
    sol, A, uA = lcpSolve(T, lcp, x, θm, stateNum, model=model, Δt=Δt)
    λ_new = sol[2stateNum+1:end]

    cons1 = model[:cons1]
    cons2 = model[:cons2]

    MOI.set.(  
        model, 
        DiffOpt.ForwardConstraintFunction(),    #to take the derivative wrt to θm
        cons1, 
        -1.0/θm^2.0*index(model[:λ][1]) + 0.0
    )

    MOI.set.(  
        model, 
        DiffOpt.ForwardConstraintFunction(),    #to take the derivative wrt to θm
        cons2, 
        [   1.0*(index(model[:v][1]) - uA[1])
            1.0*(index(model[:v][2])) - 0.0*index(model[:λ][1]) - uA[2] + lcp.sys.g*Δt     #∂cons2[2]/∂θm
        ]
    )

    if !isempty(A)
        MOI.set.(  
            model, 
            DiffOpt.ForwardObjectiveFunction(),    #to take the derivative wrt to θm
            -1.0/θm^2.0 * λ_new[1]' * A[1,1] * (ones(Float64) *index(model[:λ][1])*index(model[:λ][1])) + 0.0
        )
    end
    
    DiffOpt.forward_differentiate!(JuMP.backend(model))
    ∂sol∂θm =  MOI.get.(   
                model,
                DiffOpt.ForwardVariablePrimal(),
                [model[:q]; model[:v]; model[:λ]])

    return (sol, ∂sol∂θm)
end


function solveM(::Type{T}, lcp, x0, stateNum; Δt = 0.001, totalTimeStep = 1000) where T

    θm_actual = 0.2
    θm = 0.9
    η = 0.05

    # x0 = [0.0, rand(range(0.2, 0.5, step=0.01)), 0.1, -0.1]
    S, Λ = computeExpertTraj(lcp, x0, stateNum, θm_actual; Δt = Δt, totalTimeStep=totalTimeStep)

    for i in 1:10
        l, lg = lossGrad(S, Λ, θm, stateNum; T = T, totalTimeStep = totalTimeStep)
        println("loss = ", l, " | grad = ", lg, " | θm = ", θm)
        θm -= η*lg
    end

    return θm
end

function computeExpertTraj(lcp, x0, stateNum, θm_actual::T ; Δt = 0.001, totalTimeStep=1000) where T

    S       = Array{Array{T, 1}, 1}()
    Λ       = Array{Array{T, 1}, 1}()
    x1      = deepcopy(x0)
    push!(S, deepcopy(x1))
    #generate trajectory
    for i in 1:totalTimeStep
        expert_sol,_ ,_ = lcpSolve(T, lcp, x1, θm_actual, 2; Δt = Δt)
        x1 = expert_sol[1:2*stateNum]
        push!(S, deepcopy(x1))
        push!(Λ, deepcopy(expert_sol[2*stateNum+1:end]))
    end

    return S, Λ
end

function loss(λ, Λ; totalTimeStep = 0.001)
    l = 0.0
    for i in 1:totalTimeStep
        l += 5.0/totalTimeStep * 0.5*sum(abs2, λ[i] - Λ[i])
    end
    return l
end

function lossGrad(S, Λ, θm, stateNum; T = Float64, totalTimeStep = 1000)

    l = 0.0
    lg = 0.0

    for i in 1:totalTimeStep
        sol1, grad1 = ChainRulesCore.frule((nothing, nothing), lcpSolve, T, lcp, S[i], θm, stateNum, Δt)
        x1 = sol1[1:2*stateNum]
        λ = sol1[2*stateNum+1:end]
        l += 0.5*sum(abs2, λ - Λ[i])
        lg += (λ - Λ[i])' * grad1[2stateNum+1:end]
    end

    # loss(θm) = computeLoss(lcp, S, Λ, θm, x0, stateNum; Δt = Δt, totalTimeStep=totalTimeStep)
    # grad_m = ForwardDiff.derivative(loss, θm0)

    return l, lg
end

function computeLoss(lcp, S, Λ, θm0::T, x0, stateNum; Δt = 0.001, totalTimeStep = 1000) where T

    l = 0.0
    x1 = deepcopy(x0)
    s = Array{Array{T, 1}, 1}()
    λ = Array{Array{T, 1}, 1}()
    Q = diagm(0 => [2.0, 2.0, 0.5, 0.5])

    for i in 1:totalTimeStep

        sol,_,_ = lcpSolve(T, lcp, x1, θm0, stateNum; Δt = Δt)
        x1 = sol[1:2*stateNum]
        push!(s, deepcopy(x1))
        push!(λ, deepcopy(sol[2*stateNum+1:end]))
        l += dot(S[i] - s[i], Q*(S[i] - s[i])) + norm(Λ[i] - λ[i])
    end
    
    return l
end

function checkGradient(lcp, x1, θ1, θ2, stateNum, Δt; T = Float64)

    sol1, grad1 =  ChainRulesCore.frule((nothing, nothing), lcpSolve, T, lcp, x1, θ1, stateNum, Δt)
    q1, v1, λ1 = (sol1[1:stateNum], sol1[stateNum+1:2stateNum], sol1[2stateNum+1:end])

    sol2, grad2 =  ChainRulesCore.frule((nothing, nothing), lcpSolve, T, lcp, x1, θ2, stateNum, Δt)
    q2, v2, λ2 = (sol2[1:stateNum], sol2[stateNum+1:2stateNum], sol2[2stateNum+1:end])

    finiteDiff_q = (q2 - q1) / (θ2 - θ1)
    error_q = abs.(finiteDiff_q - grad1[1:stateNum])

    finiteDiff_v = (v2 - v1) / (θ2 - θ1)
    error_v = abs.(finiteDiff_v - grad1[stateNum+1:2stateNum])

    finiteDiff_λ = (λ2 - λ1) / (θ2 - θ1)
    error_λ = abs.(finiteDiff_λ - grad1[2stateNum+1:end])

    return error_q, error_v, error_λ, finiteDiff_q, finiteDiff_v, finiteDiff_λ, grad1
    # return error_q./finiteDiff_q, error_v./finiteDiff_v, grad1
end