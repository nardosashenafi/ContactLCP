using JuMP, Mosek, MosekTools, DiffOpt, LinearAlgebra

function lcpGradForwardControlParam(lcp, x1, θ; Δt = 0.001)

    uA      = x1[3:4]
    qA      = x1[1:2]
    qM      = qA + 0.5*Δt*uA

    x_mid   = [qM...,uA...]
    gn, γn, γt, M, h, Wn, Wt = ContactLCP.setSysAttributes(lcp, x_mid)
    ContactLCP.checkContact(lcp, gn)
    ContactLCP.trimAttributes(lcp, gn, γn, γt, M, h, Wn, Wt)

    s = lcp.current_contact_num

    if (s > 0)

        model = JuMP.Model(() -> DiffOpt.diff_optimizer(Mosek.Optimizer))
        set_silent(model)

        E = Matrix{Float64}(I, s, s)
        A = lcp.Wn'*(lcp.M\lcp.Wn)
        Minv = inv(lcp.M)

        @variable(model, λ[1:s] >= 0.0)
        @variable(model, v[1:2])
        @variable(model, q[1:2])

        b = lcp.Wn'*(lcp.M \ lcp.h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn

        @constraint(model,
            cons1,
            A*λ .+ lcp.Wn'*(Minv * ([θ * q[1], 0.0])*Δt + Minv * lcp.h*Δt) + 
            (E + diagm(0 => lcp.ϵn))*lcp.γn .>= 0.0
        )

        @constraint(model,
            cons2,
            M*(v - uA) - ([θ * q[1], 0.0])*Δt .== Wn*λ + lcp.h*Δt 
        )

        @constraint(model, 
            cons3, 
            q .== qM + 0.5*Δt*v
        )

        @objective(model, Min, dot(A*λ .+ b, λ))
        optimize!(model) 
        
        q_new = JuMP.value.(q)
        v_new = JuMP.value.(v)
        λ_new = JuMP.value.(λ)

        MOI.set.(  
            model, 
            DiffOpt.ForwardConstraintFunction(),    #to take the derivative wrt to θ
            cons2[1], 
            -q[1]*Δt         #∂cons2[1]/∂θ
        )
        
        DiffOpt.forward_differentiate!(model)
        ∂q∂θ = []
        ∂v∂θ = []

        for k in 1:2
            push!(∂q∂θ, MOI.get(   
                model,
                DiffOpt.ForwardVariablePrimal(),
                q[k]
            ))

            push!(∂v∂θ, MOI.get(   
                model,
                DiffOpt.ForwardVariablePrimal(),
                v[k]
            ))

        end

        return q_new, v_new, λ_new, ∂q∂θ, ∂v∂θ
    else
        println("Not in contact; gradient not computed for such cases yet")
    end
end

function compareFiniteDiff(lcp, func, x1, Δt, θ1, θ2)

    q1, v1, λ1, ∂q∂θ1, ∂v∂θ1 =  func(lcp, x1, θ1; Δt = Δt);
    q2, v2, λ2, _, _ =  func(lcp, x1, θ2; Δt = Δt);

    finiteDiff_q = (q2 - q1) / (θ2 - θ1)
    error_q = abs.(finiteDiff_q - ∂q∂θ1)

    finiteDiff_v = (v2 - v1) / (θ2 - θ1)
    error_v = abs.(finiteDiff_v - ∂v∂θ1)

    return error_q, error_v, error_q./finiteDiff_q, error_v./finiteDiff_v, ∂q∂θ1, ∂v∂θ1
    # return error_q./finiteDiff_q, error_v./finiteDiff_v, ∂q∂θ1, ∂v∂θ1
end

function lcpGradReverse(lcp, x1, θ; Δt = 0.001)

    uA      = x1[3:4]
    qA      = x1[1:2]
    qM      = qA + 0.5*Δt*uA

    x_mid   = [qM...,uA...]
    gn, γn, γt, M, h, Wn, Wt = ContactLCP.sysAttributes(lcp, x_mid)

    ContactLCP.checkContact(lcp, gn)
    ContactLCP.createContactMap(lcp, gn, γn, γt, M, h, Wn, Wt)

    s = lcp.current_contact_num

    if (s > 0)

        model = JuMP.Model(() -> DiffOpt.diff_optimizer(Mosek.Optimizer))
        set_silent(model)

        E = Matrix{Float64}(I, s, s)
        A = lcp.Wn'*(lcp.M\lcp.Wn)
        constraintNum = size(A, 1)

        Minv = inv(lcp.M)

        @variable(model, λ[1:s] >= 0.0)
        @variable(model, v[1:2])
        @variable(model, q[1:2])

        b = lcp.Wn'*(lcp.M \ lcp.h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn

        @constraint(model,
            cons1,
            A*λ .+ lcp.Wn'*(Minv * ([θ * q[1], 0.0])*Δt + Minv * lcp.h*Δt) + 
            (E + diagm(0 => lcp.ϵn))*lcp.γn .>= 0.0
        )

        @constraint(model,
            cons2,
            M*(v - uA) - ([θ * q[1], 0.0])*Δt .== Wn*λ + lcp.h*Δt 
        )

        @constraint(model, 
            cons3, 
            q .== qM + 0.5*Δt*v
        )

        @objective(model, Min, dot(A*λ .+ b, λ))
        optimize!(model) 
        
        q_new = JuMP.value.(q)
        v_new = JuMP.value.(v)
        λ_new = JuMP.value.(λ)

        MOI.set.(  
            model, 
            DiffOpt.ReverseVariablePrimal(),    #to take the derivative wrt to θ
            v[1]
        )
        
        DiffOpt.reverse_differentiate!(model)
        ∂q∂θ = []
        ∂v∂θ = []

        for k in 1:constraintNum
            push!(∂q∂θ, MOI.get(   
                model,
                DiffOpt.ReverseConstraintFunction(),
                cons2[constraintNum]
            ))

        end

        return q_new, v_new, λ_new, ∂q∂θ, ∂v∂θ
    else
        println("Not in contact; gradient not computed for such cases yet")
    end
end