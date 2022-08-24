mutable struct Lcp{T, TSYS}
    sys                 ::TSYS 
    total_contact_num   ::Int
    current_contact_num ::Int
    ϵn                  ::Vector{T}
    ϵt                  ::Vector{T}
    μ                   ::Vector{T}
    gn                  ::Vector{T}
    γn                  ::Vector{T}
    γt                  ::Vector{T}
    ξn                  ::Vector{T}
    ξt                  ::Vector{T}
    Wn                  ::Matrix{T}
    Wt                  ::Matrix{T}
    M                   ::Matrix{T}
    h                   ::Vector{T}

    function Lcp(sys, T)

        total_contact_num           = length(sys.contactIndex)
        current_contact_num         = sum(sys.contactIndex)
        ϵn                          = sys.ϵn
        ϵt                          = sys.ϵt
        μ                           = sys.μ
        gn                          = zeros(T, total_contact_num)
        γn                          = zeros(T, total_contact_num)
        γt                          = zeros(T, total_contact_num)
        ξn                          = zeros(T, total_contact_num)
        ξt                          = zeros(T, total_contact_num)
        gn, γn, γt, M, h, Wn, Wt    = sys(sys.x0)

        new{T, typeof(sys)}(sys, total_contact_num, current_contact_num, ϵn, ϵt, μ, gn, γn, γt, ξn, ξt, Wn, Wt, M, h)
    end
end

function checkContact(lcp::Lcp, gn)
     
    for i in 1:lcp.total_contact_num
        if gn[i] < lcp.sys.gThreshold 
            lcp.sys.contactIndex[i] = 1 
        else
            lcp.sys.contactIndex[i] = 0 
        end
    end

    lcp.current_contact_num = sum(lcp.sys.contactIndex)
end

function resetCoefficients(lcp::Lcp)
    ϵn  = lcp.sys.ϵn
    ϵt  = lcp.sys.ϵt
    μ   = lcp.sys.μ

    return ϵn, ϵt, μ
end

function sysAttributes(lcp, x)
    return lcp.sys(x)
end

function createContactMap(lcp::Lcp, gn, γn, γt, M, h, Wn, Wt)

    ϵn, ϵt, μ       = resetCoefficients(lcp)
    lcp.gn          = gn[lcp.sys.contactIndex .== 1]
    lcp.M, lcp.h    = (M, h)

    lcp.Wn          = Wn[:, lcp.sys.contactIndex .== 1]   #pick out the ones in contact
    lcp.Wt          = Wt[:, lcp.sys.contactIndex .== 1]

    #trim the coefficients and relative velocities
    lcp.ϵn, lcp.ϵt, lcp.μ, lcp.γn, lcp.γt = map(x -> x[lcp.sys.contactIndex .== 1], 
                                            [ϵn, ϵt, μ, γn, γt])
end

function solveLcp(lcp::Lcp; Δt = 0.001)

    s   = lcp.current_contact_num

    if (s > 0)
        # println("Contact detected")
        E = Matrix{Float64}(I, s, s)

        A = [lcp.Wn'*(lcp.M\(lcp.Wn - lcp.Wt*diagm(0 => lcp.μ))) lcp.Wn'*(lcp.M\lcp.Wt) zeros(s,s);
            lcp.Wt'*(lcp.M\(lcp.Wn - lcp.Wt*diagm(0 => lcp.μ))) lcp.Wt'*(lcp.M\lcp.Wt) E;
            2.0*diagm(0 => lcp.μ) -E zeros(s, s)]

        b = [lcp.Wn'*(lcp.M \ lcp.h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn;
            lcp.Wt'*(lcp.M \ lcp.h*Δt) + (E + diagm(0 => lcp.ϵt))*lcp.γt;
            zeros(s)]

        λ = lcpOpt(A, b, s)
    else
        λ = zeros(2s)
    end

    Λn = λ[1:s]
    ΛR = λ[s+1:2s]
    Λt = ΛR - diagm(0 => lcp.μ)*Λn
    
    return Λn, Λt, ΛR
end

function lcpOpt(A, b, contactNum)

    model = Model(PATHSolver.Optimizer)
    set_silent(model)
    @variable(model, λ[1:3*contactNum] >= 0.0)
    @constraints(model, begin
        (A*λ .+ b) ⟂ λ
    end)

    optimize!(model)

    return JuMP.value.(λ)
end

function oneTimeStep(lcp::Lcp, x1; Δt = 0.001)

    #TODO: replace with getstate
    uA      = x1[3:4]
    qA      = x1[1:2]
    qM      = qA + 0.5*Δt*uA

    x_mid   = [qM...,uA...]
    gn, γn, γt, M, h, Wn, Wt = sysAttributes(lcp, x_mid)

    checkContact(lcp, gn)
    createContactMap(lcp, gn, γn, γt, M, h, Wn, Wt)

    Λn, Λt, ΛR  = solveLcp(lcp; Δt=Δt)
    x2          = [qM...,uA...]

    λn = zeros(lcp.total_contact_num)
    λt = zeros(lcp.total_contact_num)
    λR = zeros(lcp.total_contact_num)
    λn[lcp.sys.contactIndex .== 1] = Λn
    λt[lcp.sys.contactIndex .== 1] = Λt
    λR[lcp.sys.contactIndex .== 1] = ΛR

    uE = lcp.M\((Wn - Wt*diagm(0 => lcp.sys.μ))*λn + Wt*λR + lcp.h*Δt) + uA
    qE = qM + 0.5*Δt*uE

    return [qE...,uE...], λn, λt

end

function fulltimestep(lcp::Lcp, T; Δt = 0.001, totalTimeStep = 1500)

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
