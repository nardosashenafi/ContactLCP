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

    function Lcp(T, sys, θ0)

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
        gn, γn, γt, M, h, Wn, Wt    = sys(sys.x0, θ0)

        new{T, typeof(sys)}(sys, total_contact_num, current_contact_num, ϵn, ϵt, μ, gn, γn, γt, ξn, ξt, Wn, Wt, M, h)
    end
end

function unstackSol(lcp::Lcp, sol)

    s = sol[1:lcp.sys.stateLength]
    λ = sol[lcp.sys.stateLength+1:end]
    return s, λ
end

#check the gap and keep track of the systems in contact
function checkContact(lcp::Lcp, gn::Vector{T}) where {T<:Real}
     
    for i in 1:lcp.total_contact_num
        if gn[i] < lcp.sys.gThreshold 
            lcp.sys.contactIndex[i] = 1 
        else
            lcp.sys.contactIndex[i] = 0 
        end
    end

    lcp.current_contact_num = sum(lcp.sys.contactIndex)
end

#the coefficients are trimmed inorder to construct A matrix and b vector. This function resets the coefficients
function resetCoefficients(lcp::Lcp)
    ϵn  = lcp.sys.ϵn
    ϵt  = lcp.sys.ϵt
    μ   = lcp.sys.μ

    return ϵn, ϵt, μ
end

function setSysAttributes(lcp, x::Vector{T}, θ::Vector{T}) where {T<:Real}
    return lcp.sys(x, θ)
end

function trimAttributes(lcp::Lcp, gn::Vector{T}, γn::Vector{T}, γt::Vector{T}, M::Matrix{T}, h::Vector{T}, Wn::Matrix{T}, Wt::Matrix{T}) where {T<:Real}

    ϵn, ϵt, μ       = resetCoefficients(lcp)
    lcp.gn          = gn[lcp.sys.contactIndex .== 1]
    lcp.M, lcp.h    = (M, h)

    lcp.Wn          = Wn[:, lcp.sys.contactIndex .== 1]   #pick out the ones in contact
    lcp.Wt          = Wt[:, lcp.sys.contactIndex .== 1]

    #trim the coefficients and relative velocities; they will be used in computing A and b 
    lcp.ϵn, lcp.ϵt, lcp.μ, lcp.γn, lcp.γt = map(x -> x[lcp.sys.contactIndex .== 1], 
                                            [ϵn, ϵt, μ, γn, γt])
end

function lcpSolve(lcp::Lcp, x::Vector{T}, θ::Vector{T}; model = JuMP.Model(Mosek.Optimizer)) where {T<:Real}

    set_silent(model)
    gn, γn, γt, M, h, Wn, Wt = setSysAttributes(lcp, x, θ)
    checkContact(lcp, gn)
    trimAttributes(lcp, gn, γn, γt, M, h, Wn, Wt)

    E = Matrix{T}(I, lcp.current_contact_num, lcp.current_contact_num)

    Minv = inv(lcp.M)
    A = lcp.Wn'*(Minv*lcp.Wn)
    b = lcp.Wn'*(Minv*lcp.h*lcp.sys.Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn

    qM, uA = lcp.sys(x)

    @variable(model, λ[1:lcp.total_contact_num] >= T(0.0))
    @variable(model, v[1:length(uA)])
    @variable(model, q[1:length(qM)])

    @expression(model, contactForces, [lcp.sys.contactIndex[i] == 1 ? λ[i] : T(0.0) for i in 1:lcp.total_contact_num] )
    @constraint(model,
        cons1,
        A*λ[lcp.sys.contactIndex .== 1] .+ b .>= T(0.0)
    )

    @constraint(model,
        cons2,
        M*(v - uA) - h*lcp.sys.Δt .== Wn*contactForces
    )

    @constraint(model, 
        cons3, 
        q .== qM + 0.5*lcp.sys.Δt*v
    )

    @objective(model, Min, dot(A*λ[lcp.sys.contactIndex .== 1] .+ b, λ[lcp.sys.contactIndex .== 1]))

    optimize!(model) 

   return vcat(JuMP.value.(q), JuMP.value.(v), JuMP.value.(λ)), A, b
end

function oneTimeStep(lcp::Lcp, x1::Vector{T}, θ::Vector{T}) where {T<:Real}

    qA, uA              = lcp.sys(x1)
    qM                  = qA + 0.5*lcp.sys.Δt*uA
    
    x_mid               = [qM...,uA...]
    sol, _, _           = lcpSolve(lcp, x_mid, θ)
    s, λn               = unstackSol(lcp, sol)
    qE, uE              = lcp.sys(s)
    return [qE...,uE...], λn

end

function fulltimestep(lcp::Lcp, x0::Vector{T}, θ::Vector{T}) where {T<:Real}

    X       = Vector{Vector{Float64}}()
    Λn      = Vector{Vector{Float64}}()
    t       = Vector{T}()

    if isempty(x0)
        x = deepcopy(lcp.sys.x0)
    else
        x = deepcopy(x0)
    end

    X       = push!(X, x)
    t       = push!(t, T(0.0))

    for i in 1:lcp.sys.totalTimeStep
        x, λn = oneTimeStep(lcp, x, θ)
        push!(X, x)
        push!(Λn, λn)
        push!(t, t[end]+lcp.sys.Δt)
    end

    return X, t, Λn
end
