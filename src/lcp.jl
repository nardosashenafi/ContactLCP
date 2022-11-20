struct Lcp{T, TSYS}
    sys ::TSYS 
    function Lcp(T, sys)
        new{T, typeof(sys)}(sys)
    end
end

function sysAttributes(lcp::Lcp, x, param; kwargs...)
    return lcp.sys(x, param; kwargs...)
end

function checkContact(gn::Vector{T}, gThreshold, total_contact_num) where {T<:Real}
     
    contactIndex = zeros(T, total_contact_num)

    for i in 1:total_contact_num
        if (gn[i] < gThreshold) 
            contactIndex[i] = 1 
        end
    end
    current_contact_num = Int(sum(contactIndex))

    return contactIndex, current_contact_num
end

function trimAttributes(gn, γn, γt, Wn, Wt, ϵn, ϵt, μ, contactIndex)

    gnt         = gn[contactIndex .== 1]
    Wnt         = Wn[:, contactIndex .== 1]   #pick out the ones in contact
    Wtt         = Wt[:, contactIndex .== 1]

    #trim the coefficients and relative velocities
    ϵnt, ϵtt, μt, γnt, γtt = map(x -> x[contactIndex .== 1], 
                                [ϵn, ϵt, μ, γn, γt])

    return gnt, ϵnt, ϵtt, μt, γnt, γtt, Wnt, Wtt
end

function getAb(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, contactIndex, current_contact_num; Δt = 0.001f0)

    gnt, ϵnt, ϵtt, μt, γnt, γtt, Wnt, Wtt = trimAttributes(gn, γn, γt, Wn, Wt, ϵn, ϵt, μ, contactIndex)
    s = current_contact_num
    E = Matrix{Float32}(I, s, s)

    A = [Wnt'*(M\(Wnt - Wtt*diagm(0 => μt))) Wnt'*(M\Wtt) zeros(Float32, s,s);
        Wtt'*(M\(Wnt - Wtt*diagm(0 => μt))) Wtt'*(M\Wtt) E;
        2.0f0*diagm(0 => μt) -E zeros(Float32, s, s)]

    b = [Wnt'*(M \ h*Δt) + (E + diagm(0 => ϵnt))*γnt;
        Wtt'*(M \ h*Δt) + (E + diagm(0 => ϵtt))*γtt;
        zeros(Float32, s)]

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

function solveLcp(gn, γn, γt, M, h::Vector{T}, Wn, Wt, ϵn, ϵt, μ, gThreshold; Δt=0.001f0) where {T<:Real}

    total_contact_num = length(gn)
    contactIndex, s = checkContact(gn, gThreshold, total_contact_num)

    Λn           = zeros(T, total_contact_num)
    ΛR           = zeros(T, total_contact_num)
    Λt           = zeros(T, total_contact_num)

    if s > 0
        A, b = getAb(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, contactIndex, s; Δt = Δt)
        λ    =  lemkeLexi(A, b, [])
        # λ = lcpOpt(A, b, s)

        Λn[contactIndex .== 1] = λ[1:s]
        ΛR[contactIndex .== 1] = λ[s+1:2s]
        Λt[contactIndex .== 1] = λ[s+1:2s] - diagm(0 => μ[contactIndex .== 1])*λ[1:s]
    end

    return Λn, Λt, ΛR
end


function solveLcp(gn, γn, γt, M, h::Vector{T}, Wn, Wt, ϵn, ϵt, μ, gThreshold, x; Δt=0.001f0) where {T<:Real}

    total_contact_num = length(gn)
    contactIndex, s = checkContact(gn, gThreshold, total_contact_num)

    Λn           = zeros(T, total_contact_num)
    ΛR           = zeros(T, total_contact_num)
    Λt           = zeros(T, total_contact_num)

    if s > 0
        A, b = getAb(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, contactIndex, s; Δt = Δt)
        λ =  lemkeLexi(A, b, x)
        # λ = lcpOpt(A, b, s)

        Λn[contactIndex .== 1] = λ[1:s]
        ΛR[contactIndex .== 1] = λ[s+1:2s]
        Λt[contactIndex .== 1] = λ[s+1:2s] - diagm(0 => μ[contactIndex .== 1])*λ[1:s]
    end

    return Λn, Λt, ΛR
end

function oneTimeStep(lcp::Lcp, x1, param::Vector{T}; Δt = 0.001, kwargs...) where {T<:Real}

    qA, uA  = lcp.sys(x1)
    qM      = qA + 0.5f0*Δt*uA

    x_mid   = vcat(qM, uA)
    gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold = sysAttributes(lcp, x_mid, param; kwargs...)
    λn, λt, λR  = solveLcp(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold, x1; Δt=Δt)

    uE = M\((Wn - Wt*diagm(0 => μ))*λn + Wt*λR + h*Δt) + uA
    qE = qM + 0.5f0*Δt*uE

    return vcat(qE,uE), λn, λt
end

function fulltimestep(sys, x0, param::Vector{T}, total_contact_num; Δt = 0.001f0, totalTimeStep = 500) where {T<:Real}
    lcp = Lcp(T, sys)
    fulltimestep(lcp, x0, param, total_contact_num; Δt = Δt, totalTimeStep = totalTimeStep) 
end

function fulltimestep(lcp::Lcp, x0, param::Vector{T}; Δt = 0.001f0, totalTimeStep = 500, kwargs...) where {T<:Real}

    X       = Vector{Vector{T}}(undef, totalTimeStep)
    Λn      = Vector{Vector{T}}(undef, totalTimeStep)
    Λt      = Vector{Vector{T}}(undef, totalTimeStep)
    t       = Vector{T}(undef, totalTimeStep)
    x       = deepcopy(x0)
    ti      = 0.0f0

    oneStep(x, param) = oneTimeStep(lcp, x, param; Δt=Δt, kwargs...)

    for i in 1:totalTimeStep
        x, λn, λt  = oneStep(x, param)
        X[i]    = x
        Λn[i]   = λn
        Λt[i]   = λt
        ti      = ti+Δt
        t[i]    = ti
    end

    return X, t, Λn, Λt
end
