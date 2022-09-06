mutable struct ContactMap{T, TSYS}
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

    function ContactMap(T, sys, θ0)

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

function unstackSol(cm::ContactMap, sol)

    s = sol[1:cm.sys.stateLength]
    λ = sol[cm.sys.stateLength+1:end]
    return s, λ
end

#check the gap and keep track of the systems in contact
function checkContact(cm::ContactMap, gn::Vector{T}) where {T<:Real}
     
    for i in 1:cm.total_contact_num
        if gn[i] < cm.sys.gThreshold 
            cm.sys.contactIndex[i] = 1 
        else
            cm.sys.contactIndex[i] = 0 
        end
    end

    cm.current_contact_num = sum(cm.sys.contactIndex)
end

function setSysAttributes(cm, x::Vector{T}, θ::Vector{T}) where {T<:Real}
    return cm.sys(x, θ)
end

function cmSolve(cm::ContactMap, x::Vector{T}, θ::Vector{T}; model = JuMP.Model(Mosek.Optimizer)) where {T<:Real}

    set_silent(model)
    gn, γn, γt, M, h, Wn, Wt = setSysAttributes(cm, x, θ)
    checkContact(cm, gn)
    trimAttributes(cm, gn, γn, γt, M, h, Wn, Wt)

    E = Matrix{T}(I, cm.current_contact_num, cm.current_contact_num)

    Minv = inv(cm.M)
    A = cm.Wn'*(Minv*cm.Wn)
    b = cm.Wn'*(Minv*cm.h*cm.sys.Δt) + (E + diagm(0 => cm.ϵn))*cm.γn

    qM, uA = cm.sys(x)

    @variable(model, λ[1:cm.total_contact_num] >= T(0.0))
    @variable(model, v[1:length(uA)])
    @variable(model, q[1:length(qM)])

    @expression(model, contactForces, [cm.sys.contactIndex[i] == 1 ? λ[i] : T(0.0) for i in 1:cm.total_contact_num] )
    @constraint(model,
        cons1,
        A*λ[cm.sys.contactIndex .== 1] .+ b .>= T(0.0)
    )

    @constraint(model,
        cons2,
        M*(v - uA) - h*cm.sys.Δt .== Wn*contactForces
    )

    @constraint(model, 
        cons3, 
        q .== qM + 0.5*cm.sys.Δt*v
    )

    @objective(model, Min, dot(A*λ[cm.sys.contactIndex .== 1] .+ b, λ[cm.sys.contactIndex .== 1]))

    optimize!(model) 

   return vcat(JuMP.value.(q), JuMP.value.(v), JuMP.value.(λ)), A, b
end

function oneTimeStep(cm::ContactMap, x1::Vector{T}, θ::Vector{T}) where {T<:Real}

    uA                  = x1[3:4]
    qA                  = x1[1:2]
    qM                  = qA + 0.5*cm.sys.Δt*uA
    
    x_mid               = [qM...,uA...]

    gn, γn, γt, M, h, Wn, Wt = setSysAttributes(cm, x_mid, θ)
    checkContact(cm, gn)

    if any(cm.sys.contactIndex .== 1)

        Ξ = impactMap(sys, ϕ)
        uE = Ξ*uA
        #switch θ to the next spoke
        qM[1] = -qM[1]
        qE .== qM + 0.5*cm.sys.Δt*uE

    else

        uE = inv(M)*(h*cm.sys.Δt) + uA
        qE .== qM + 0.5*cm.sys.Δt*uE

    end
    
    return [qE...,uE...]

end

function fulltimestep(cm::ContactMap, x0::Vector{T}, θ::Vector{T}) where {T<:Real}

    X       = Vector{Vector{Float64}}()
    Λn      = Vector{Vector{Float64}}()
    t       = Vector{T}()

    if isempty(x0)
        x = deepcopy(cm.sys.x0)
    else
        x = deepcopy(x0)
    end

    X       = push!(X, x)
    t       = push!(t, T(0.0))

    for i in 1:cm.sys.totalTimeStep
        x, λn = oneTimeStep(cm, x, θ)
        push!(X, x)
        push!(Λn, λn)
        push!(t, t[end]+cm.sys.Δt)
    end

    return X, t, Λn
end
