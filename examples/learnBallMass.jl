
function lcpGradForwardSysParam(lcp, x1, θm; Δt = 0.001)

    uA      = x1[3:4]
    qA      = x1[1:2]
    qM      = qA + 0.5*Δt*uA

    x_mid   = [qM...,uA...]
    gn, γn, γt, M, h, Wn, Wt = ContactLCP.sysAttributes(lcp, x_mid)

    ContactLCP.checkContact(lcp, gn)
    ContactLCP.createContactMap(lcp, gn, γn, γt, M, h, Wn, Wt)

    s = lcp.current_contact_num

    #change the mass, update the mass matrix and the gravitational force
    M = [θm 0.0; 0.0 θm]
    Minv = inv(M)
    h = [0.0, -θm*sys.g]

    if (s > 0)

        model = JuMP.Model(() -> DiffOpt.diff_optimizer(Mosek.Optimizer))
        set_silent(model)

        E = Matrix{Float64}(I, s, s)
        A = lcp.Wn'*(M\lcp.Wn)

        @variable(model, λ[1:s] >= 0.0)
        @variable(model, v[1:2])
        @variable(model, q[1:2])

        b = lcp.Wn'*(M \ h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn

        @constraint(model,
            cons1,
            A*λ .+ lcp.Wn'*(Minv * h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn .>= 0.0
        )

        @constraint(model,
            cons2,
            (M[i,i]*(v[i] - uA[i]) for i in 1:2) .== Wn*λ + h*Δt 
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
            DiffOpt.ForwardConstraintFunction(),    #to take the derivative wrt to θm
            cons2[1], 
            (v[1] - uA[1])      #∂cons2[1]/∂θm
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