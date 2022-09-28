mutable struct Lcp{T, TSYS}
    sys                 ::TSYS 
    total_contact_num   ::Int
    current_contact_num ::Int
    ϵn                  ::Vector{T}
    ϵt                  ::Vector{T}
    μ                   ::Vector{T}

    function Lcp(sys, T)

        total_contact_num           = length(sys.contactIndex)
        current_contact_num         = sum(sys.contactIndex)
        ϵn                          = sys.ϵn
        ϵt                          = sys.ϵt
        μ                           = sys.μ

        new{T, typeof(sys)}(sys, total_contact_num, current_contact_num, ϵn, ϵt, μ)
    end
end

function sysAttributes(lcp, x)
    return lcp.sys(x, [])
end

function checkContact(lcp::Lcp, gn, γn)
     
    contactIndex = zeros(lcp.total_contact_num)

    for i in 1:lcp.total_contact_num
        if (gn[i] < lcp.sys.gThreshold) && (γn[i] < 0.0)
            contactIndex[i] = 1 
        end
    end

    lcp.current_contact_num = sum(contactIndex)
    return contactIndex
end

function resetCoefficients(lcp::Lcp)
    return lcp.sys.ϵn, lcp.sys.ϵt, lcp.sys.μ
end

function trimAttributes(lcp::Lcp, gn, γn, γt, Wn, Wt)

    contactIndex = checkContact(lcp, gn, γn)

    ϵn, ϵt, μ   = resetCoefficients(lcp)
    gnt         = gn[contactIndex .== 1]
    Wnt         = Wn[:, contactIndex .== 1]   #pick out the ones in contact
    Wtt         = Wt[:, contactIndex .== 1]

    #trim the coefficients and relative velocities
    ϵnt, ϵtt, μt, γnt, γtt = map(x -> x[contactIndex .== 1], 
                                [ϵn, ϵt, μ, γn, γt])

    return gnt, ϵnt, ϵtt, μt, γnt, γtt, Wnt, Wtt
end

function getAb(lcp::Lcp, gn, γn, γt, M, h, Wn, Wt; Δt = 0.001)

    gnt, ϵnt, ϵtt, μt, γnt, γtt, Wnt, Wtt = trimAttributes(lcp, gn, γn, γt, Wn, Wt)

    s   = lcp.current_contact_num

    # println("Contact detected")
    E = Matrix{Float64}(I, s, s)

    A = [Wnt'*(M\(Wnt - Wtt*diagm(0 => μt))) Wnt'*(M\Wtt) zeros(s,s);
        Wtt'*(M\(Wnt - Wtt*diagm(0 => μt))) Wtt'*(M\Wtt) E;
        2.0*diagm(0 => μt) -E zeros(s, s)]

    b = [Wnt'*(M \ h*Δt) + (E + diagm(0 => ϵnt))*γnt;
        Wtt'*(M \ h*Δt) + (E + diagm(0 => ϵtt))*γtt;
        zeros(s)]

    return A, b
end

function lcpOpt(A, b, contactNum)

    if contactNum > 0
        model = Model(PATHSolver.Optimizer)
        set_silent(model)
        @variable(model, λ[1:3*contactNum] >= 0.0)
        @constraints(model, begin
            (A*λ .+ b) ⟂ λ
        end)

        optimize!(model)

        return JuMP.value.(λ)
    else 
        return zeros(contactNum)
    end
end

function solveLcp(lcp, gn, γn, γt, M, h, Wn, Wt; Δt=0.001)

    contactIndex = checkContact(lcp, gn, γn)
    s            = lcp.current_contact_num
    λn = zeros(lcp.total_contact_num)
    λt = zeros(lcp.total_contact_num)
    λR = zeros(lcp.total_contact_num)

    if s > 0
        A, b         = getAb(lcp, gn, γn, γt, M, h, Wn, Wt; Δt = Δt)
        # λ  = lcpOpt(A, b, s)
        λ =  lemke(A, b)

        Λn = λ[1:s]
        ΛR = λ[s+1:2s]
        Λt = ΛR - diagm(0 => lcp.sys.μ[contactIndex .== 1])*Λn

        λn[contactIndex .== 1] = Λn
        λt[contactIndex .== 1] = Λt
        λR[contactIndex .== 1] = ΛR
    end

    return λn, λt, λR
end

function solveLcp(lcp::Lcp, x; Δt = 0.001)
    gn, γn, γt, M, h, Wn, Wt = sysAttributes(lcp, x)
    return solveLcp(lcp, gn, γn, γt, M, h, Wn, Wt; Δt=0.001)
end

function oneTimeStep(lcp::Lcp, x1; Δt = 0.001)

    #TODO: replace with getstate
    qA, uA  = lcp.sys(x1)
    qM      = qA + 0.5f0*Δt*uA

    x_mid   = vcat(qM, uA)
    gn, γn, γt, M, h, Wn, Wt = sysAttributes(lcp, x_mid)
    # println("x = ", x_mid)
    λn, λt, λR  = solveLcp(lcp, gn, γn, γt, M, h, Wn, Wt; Δt=Δt)
    x2          = vcat(qM,uA)

    uE = M\((Wn - Wt*diagm(0 => lcp.sys.μ))*λn + Wt*λR + h*Δt) + uA
    qE = qM + 0.5f0*Δt*uE

    return vcat(qE,uE), λn, λt

end

function fulltimestep(lcp::Lcp, T; Δt = 0.001, totalTimeStep = 500)

    X       = Array{Array{T, 1}, 1}()
    Λn      = Array{Array{T, 1}, 1}()
    Λt      = Array{Array{T, 1}, 1}()
    t       = Array{T, 1}()
    x       = deepcopy(lcp.sys.x0)
    X       = push!(X, x)
    t       = push!(t, 0.0)

    for i in 1:totalTimeStep
        x, λn, λt  = oneTimeStep(lcp, x; Δt=Δt)
        push!(X, x)
        push!(Λn, λn)
        push!(Λt, λt)
        push!(t, t[end]+Δt)
    end

    return X, t, Λn, Λt
end
