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

function setSysAttributes(cm, x, θ::Vector{T}; limitcycle=false) where {T<:Real}
    return cm.sys(x, θ; limitcycle=limitcycle)
end

function totalEnergy(sys, x)
    θ, ϕ = x[1:2]
    θdot, ϕdot = x[3:4]
    return 0.5*(sys.I1 + sys.mt*sys.l1^2)*θdot^2 - sys.m2*sys.l1*sys.l2*θdot*ϕdot*cos(θ - ϕ) + 
            0.5*(sys.I2+sys.m2*sys.l2^2)*ϕdot^2 + sys.mt*sys.g*sys.l1*cos(θ - sys.γ) - sys.m2*sys.g*sys.l2*cos(ϕ - sys.γ)
end

function limitCycle(cm::ContactMap) 

    #limitcycle with no slope
    γ               = deepcopy(cm.sys.γ)
    cm.sys.γ        = 0.0
    totalTimeStep   = deepcopy(cm.sys.totalTimeStep)
    cm.sys.totalTimeStep = 5000

    x0      = [0.2, 0.1, -2.0, 0.0]
    X, t    = fulltimestep(cm, x0, [100.0, 20.0]; limitcycle=true)

    #reset γ back for the rest of the computations
    cm.sys.γ = deepcopy(γ)
    cm.sys.totalTimeStep = deepcopy(totalTimeStep)

    x = X[end-1]
    t = t[end-1] .- t[end-1][1]

    figure()
    plot(getindex.(x, 1), getindex.(x, 3))
    ylabel("Limit Cycle check", fontsize=15)
    
    return x, t

end

function fulltimestep(cm::ContactMap, x1, θ::Vector{T}; limitcycle=false, timeSteps = cm.sys.totalTimeStep, X = Vector{Vector{Vector{T}}}(), t = Vector{Vector{T}}()) where {T<:Real}

    if isempty(x0)
        x = deepcopy(cm.sys.x0)
    else
        x = deepcopy(x0)
    end
    push!(X, deepcopy([x1]))
    push!(t, T.([0.0]))

    for i in 1:timeSteps
        qA, uA      = cm.sys(x1)
        qM          = qA + 0.5*cm.sys.Δt*uA
        
        x_mid       = [qM...,uA...]

        gn, γn, γt, M, h, Wn, Wt = setSysAttributes(cm, x_mid, θ; limitcycle=limitcycle)
        checkContact(cm, gn)

        if any(cm.sys.contactIndex .== 1) && any(vnormal(cm.sys, x_mid) .< 0.0)

            ϕ       = x_mid[2]
            Ξ       = impactMap(sys, ϕ)
            uE      = Ξ*uA
            # println("Preimpact KE = ", uA'*cm.M*uA, " postimpact KE = ", uE'*cm.M*uE)

            #switch θ to the next spoke
            if uA[1] < 0.0
                qM[1] = cm.sys.α
            elseif uA[1] > 0.0
                qM[1] = -cm.sys.α
            end
            qE = qM + 0.5*cm.sys.Δt*uE
            x2 = [qE...,uE...]
            push!(X, deepcopy([x2]))
            push!(t, [t[end][end]+cm.sys.Δt])

            # println("Energy preimpact = ", totalEnergy(cm.sys, x1), " postimpact = ", totalEnergy(cm.sys, [qE...,uE...]))
        else

            uE = inv(M)*(h*cm.sys.Δt) + uA
            qE = qM + 0.5*cm.sys.Δt*uE
            x2 = [qE...,uE...]
            push!(X[end], deepcopy(x2))
            push!(t[end], t[end][end]+cm.sys.Δt)
        end
        x1 = x2
    end

    return X, t
end

