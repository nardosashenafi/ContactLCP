mutable struct Lcp{T, TSYS}
    sys                 ::TSYS 
    total_contact_num   ::Int
    current_contact_num ::Int
    ϵn                  ::Vector{T}
    ϵt                  ::Vector{T}
    μ                   ::Vector{T}

    function Lcp(T, sys)

        total_contact_num           = length(sys.contactIndex)
        current_contact_num         = sum(sys.contactIndex)
        ϵn                          = sys.ϵn
        ϵt                          = sys.ϵt
        μ                           = sys.μ

        new{T, typeof(sys)}(sys, total_contact_num, current_contact_num, ϵn, ϵt, μ)
    end
end

function sysAttributes(lcp::Lcp, x, param)
    return lcp.sys(x, param)
end

function checkContact(lcp::Lcp, gn, γn)
     
    contactIndex = zeros(lcp.total_contact_num)

    for i in 1:lcp.total_contact_num
        # if (gn[i] < lcp.sys.gThreshold) && (γn[i] <= 0.0)
        if (gn[i] < lcp.sys.gThreshold) 
            contactIndex[i] = 1 
        end
    end

    lcp.current_contact_num = sum(contactIndex)
    return contactIndex
end

#the coefficients are trimmed inorder to construct A matrix and b vector. This function resets the coefficients
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

function solveLcp(lcp::Lcp, gn::Vector{T}, γn, γt, M, h, Wn, Wt; Δt=0.001) where {T<:Real}

    contactIndex = checkContact(lcp, gn, γn)
    # println("Gn = ", gn[4:5])
    # println("γn = ", γn[4:5])
    println("ContactIndex = ", contactIndex)
    s            = lcp.current_contact_num
    Λn           = zeros(T, lcp.total_contact_num)
    ΛR           = zeros(T, lcp.total_contact_num)
    Λt           = zeros(T, lcp.total_contact_num)

    if s > 0
        A, b = getAb(lcp, gn, γn, γt, M, h, Wn, Wt; Δt = Δt)
        # λ  = lcpOpt(A, b, s)
        λ =  lemke(A, b)

        Λn[contactIndex .== 1] = λ[1:s]
        ΛR[contactIndex .== 1] = λ[s+1:2s]
        Λt[contactIndex .== 1] = λ[s+1:2s] - diagm(0 => lcp.sys.μ[contactIndex .== 1])*λ[1:s]
    end

    return Λn, Λt, ΛR
end

function solveLcp(lcp::Lcp, x, param; Δt = 0.001)
    gn, γn, γt, M, h, Wn, Wt = sysAttributes(lcp, x, param)
    return solveLcp(lcp, gn, γn, γt, M, h, Wn, Wt; Δt=0.001)
end

function oneTimeStep(lcp::Lcp, x1, param; Δt = 0.001)

    qA, uA  = lcp.sys(x1)
    qM      = qA + 0.5f0*Δt*uA

    x_mid   = vcat(qM, uA)
    gn, γn, γt, M, h, Wn, Wt = sysAttributes(lcp, x_mid, param)
    println("x = ", x_mid)
    λn, λt, λR  = solveLcp(lcp, gn, γn, γt, M, h, Wn, Wt; Δt=Δt)
    x2          = vcat(qM,uA)

    uE = M\((Wn - Wt*diagm(0 => lcp.sys.μ))*λn + Wt*λR + h*Δt) + uA
    qE = qM + 0.5f0*Δt*uE

    return vcat(qE,uE), λn, λt
end

function fulltimestep(lcp::Lcp, x0::Vector{T}, param; Δt = 0.001, totalTimeStep = 500) where {T<:Real}

    X       = Vector{Vector{T}}(undef, totalTimeStep+1)
    Λn      = Vector{Vector{T}}(undef, totalTimeStep+1)
    Λt      = Vector{Vector{T}}(undef, totalTimeStep+1)
    t       = Vector{T}(undef, totalTimeStep+1)
    x       = deepcopy(x0)
    X[1]    = x
    t[1]    = 0.0
    Λn[1]   = zeros(T, lcp.total_contact_num)
    Λt[1]   = zeros(T, lcp.total_contact_num)

    for i in 2:totalTimeStep+1
        x, λn, λt  = oneTimeStep(lcp, x, param; Δt=Δt)
        X[i]    = x
        Λn[i]   = λn
        Λt[i]   = λt
        t[i]    = t[i-1]+Δt
    end

    return X, t, Λn, Λt
end
