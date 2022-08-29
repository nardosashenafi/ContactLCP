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

function setSysAttributes(lcp, x)
    return lcp.sys(x)
end

function trimAttributes(lcp::Lcp, gn, γn, γt, M, h, Wn, Wt)

    ϵn, ϵt, μ       = resetCoefficients(lcp)
    lcp.gn          = gn[lcp.sys.contactIndex .== 1]
    lcp.M, lcp.h    = (M, h)

    lcp.Wn          = Wn[:, lcp.sys.contactIndex .== 1]   #pick out the ones in contact
    lcp.Wt          = Wt[:, lcp.sys.contactIndex .== 1]

    #trim the coefficients and relative velocities; they will be used in computing A and b 
    lcp.ϵn, lcp.ϵt, lcp.μ, lcp.γn, lcp.γt = map(x -> x[lcp.sys.contactIndex .== 1], 
                                            [ϵn, ϵt, μ, γn, γt])
end

function solveLcp(lcp::Lcp, x; Δt = 0.001)
    return lcpOptNormal(lcp, x, 2; Δt = Δt)
end

function lcpOptNormal(lcp::Lcp, x, stateNum; Δt = 0.001)   
    #this only considers normal contact forces. For tangential forces, use PATHSolver

    model = Model(Mosek.Optimizer)
    set_silent(model)

    gn, γn, γt, M, h, Wn, Wt = setSysAttributes(lcp, x)
    checkContact(lcp, gn)
    trimAttributes(lcp, gn, γn, γt, M, h, Wn, Wt)

    E = Matrix{Float64}(I, lcp.current_contact_num, lcp.current_contact_num)

    A = lcp.Wn'*(lcp.M\lcp.Wn)
    b = lcp.Wn'*(lcp.M \ lcp.h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn

    #TODO: use env macro to directly receive state info from the system struct
    qM = x[1:stateNum]
    uA = x[stateNum+1:2stateNum]

    @variable(model, λ[1:lcp.total_contact_num] >= 0.0)
    @variable(model, v[1:stateNum])
    @variable(model, q[1:stateNum])

    @expression(model, contactForces, [lcp.sys.contactIndex[i] == 1 ? λ[i] : 0.0 for i in 1:lcp.total_contact_num] )
    @constraints(model, begin
        A*λ[lcp.sys.contactIndex .== 1] .+ b .>= 0.0
        M*(v - uA) .== Wn*contactForces + h*Δt
        q .== qM + 0.5*Δt*v
    end)

    @objective(model, Min, dot(A*λ[lcp.sys.contactIndex .== 1] .+ b, λ[lcp.sys.contactIndex .== 1]))
    optimize!(model)

    return JuMP.value.(q), JuMP.value.(v), JuMP.value.(λ)
end

function oneTimeStep(lcp::Lcp, x1; Δt = 0.001)

    #TODO: replace with getstate
    uA      = x1[3:4]
    qA      = x1[1:2]
    qM      = qA + 0.5*Δt*uA

    x_mid   = [qM...,uA...]
    qE, uE, λn  = solveLcp(lcp, x_mid; Δt=Δt)

    return [qE...,uE...], λn

end

function fulltimestep(lcp::Lcp, T; Δt = 0.001, totalTimeStep = 1500)

    X       = Array{Array{T, 1}, 1}()
    Λn      = Array{Array{T, 1}, 1}()
    t       = Array{T, 1}()
    x       = deepcopy(lcp.sys.x0)
    X       = push!(X, x)
    t       = push!(t, 0.0)

    for i in 1:totalTimeStep
        x, λn = oneTimeStep(lcp, x; Δt=Δt)
        push!(X, x)
        push!(Λn, λn)
        push!(t, t[end]+Δt)
    end

    return X, t, Λn
end
