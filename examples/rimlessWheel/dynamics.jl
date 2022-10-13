export RimlessWheel

mutable struct RimlessWheel{T}  
    m1              ::T         #wheel 
    m2              ::T         #torso
    I1              ::T 
    I2              ::T
    mt              ::T
    l1              ::T
    l2              ::T 
    g               ::T 
    k               ::Int
    α               ::T
    γ               ::T
    ϵn              ::Vector{T}
    ϵt              ::Vector{T}
    μ               ::Vector{T}
    x0              ::Vector{T}
    contactIndex    ::Vector{T}
    gThreshold      ::T

    function RimlessWheel(T)

        m1              = 2.32     
        m2              = 4.194 
        I1              = 0.0784160
        I2              = 0.0380256
        mt              = m1 + m2
        l1              = 0.26          #wheel
        l2              = 0.05          #torso
        g               = T(9.81)
        k               = 10
        α               = T(360.0/k/2.0 * pi/180.0)
        γ               = T(0.0*pi/180.0)
        ϵn              = T.(0.0*ones(k))
        ϵt              = T.(0.0*ones(k))
        μ               = T.(0.6*ones(k))
        x0              = zeros(T, 8)
        #[x, y, ϕ, θ, xdot, ydot, ϕdot, θdot]
        # state x is parallel to the plane, y points normal to the plane
        # θ is that of the spoke in contact with the plane initially
        contactIndex    = zeros(T, k)
        gThreshold      = T.(0.001)

        new{T}(m1, m2, I1, I2, mt, l1, l2, g, k, α, γ, ϵn, ϵt, μ, x0, contactIndex, gThreshold)
    end

end

function initialState(sys, θ0, θ0dot, ϕ0, ϕ0dot)
    @assert -sys.α <= θ0 <= sys.α "Give an initial spoke angle for the spoke in contact"
    x = sys.l1*sin(θ0)
    y = sys.l1*cos(θ0)
    ϕ = ϕ0 
    θ = θ0 
    xdot = -sys.l1*cos(θ0) * θ0dot
    ydot = sys.l1*sin(θ0) * θ0dot
    ϕdot = ϕ0dot 
    θdot = θ0dot

    return [x, y, ϕ, θ, xdot, ydot, ϕdot, θdot]
end

#returns the attributes need to model contact
function (sys::RimlessWheel)(x, θ::Vector{T}; expert=false) where {T<:Real}
    gn  = gap(sys, x)  
    γn  = vnormal(sys, x)
    γt  = vtang(sys, x)
    M   = massMatrix(sys, x)
    h   = genForces(sys, x, θ; expert=expert)
    Wn  = wn(sys, x)
    Wt  = wt(sys, x)

    return gn, γn, γt, M, h, Wn, Wt
end

#parses out the poses and velocities from the state vector
function (sys::RimlessWheel)(x::Vector{T}) where {T<:Real}
    q = x[1:4]      #[x, y, ϕ, θ]
    u = x[5:8]      #[xdot, ydot, ϕdot, θdot]
    return q, u
end

function setCoefficients(sys::RimlessWheel, ϵn, ϵt, μ)
    sys.ϵn  = ϵn
    sys.ϵt  = ϵt
    sys.μ   = μ
end

function setInitial(sys::RimlessWheel, z)
    sys.x0 = x
end

function gap(sys::RimlessWheel, z)
    
    k_range = range(1, stop=sys.k, step=1)
    y, θ = (z[2], z[4])
    return y .+ sys.l1*cos.(θ .+ 2*sys.α*k_range) 
end

function vnormal(sys::RimlessWheel, z)
    Wn = wn(sys, z)
    _, v = sys(z)
    return Wn'*v
end

function vtang(sys::RimlessWheel, z)
    Wt = wt(sys, z)
    _, v = sys(z)
    return Wt'*v
end

function massMatrix(sys::RimlessWheel, z)
    q, v = sys(z)
    x, y, ϕ, θ = q
    
    M = [sys.mt 0.0 sys.m2*sys.l2*cos(ϕ) 0.0;
        0.0 sys.mt sys.m2*sys.l2*sin(ϕ) 0.0;
        sys.m2*sys.l2*cos(ϕ) sys.m2*sys.l2*sin(ϕ) sys.I2+sys.m2*sys.l2^2 0.0;
        0.0 0.0 0.0 sys.I1]

    return M
end

inputLayer(x) = [cos(x[1]), sin(x[1]), cos(x[2]), sin(x[2]), x[3], x[4]]

function control(z, θp::Vector{T}; expert=false, satu = Inf) where {T<:Real}
    q, v = sys(z)

    if expert #working expert controller
        return -θp[1]*(q[3]-0.38f0) - θp[2]*v[3]
    else
        u = -θp[1]*(q[3]-0.38f0) - θp[2]*v[3]
        return clamp(u, -satu, satu)
        # return 0.0
        # @assert length(θp) == 6 + DiffEqFlux.paramlength(Hd)
        # y = MLBasedESC.controller(npbc, inputLayer(x), θp)
        # y = unn(inputLayer(x), θp)[1]
        # return clamp(y, -satu, satu)
    end
end

function genForces(sys::RimlessWheel, z, param::Vector{T}; expert=false) where {T<:Real}

    q, v = sys(z)
    x, y, ϕ, θ = q
    xdot, ydot, ϕdot, θdot = v
    #h = Bu - C qdot - G
    B = [0.0f0, 0.0f0, 1.0f0, -1.0f0]

    C = [0.0f0 0.0f0 -sys.m2*sys.l2*sin(ϕ)*ϕdot 0.0f0;
        0.0f0 0.0f0 sys.m2*sys.l2*cos(ϕ)*ϕdot 0.0f0;
        -sys.m2*sys.l2*sin(ϕ)*ϕdot sys.m2*sys.l2*cos(ϕ)*ϕdot 0.0f0 0.0f0;
        0.0f0 0.0f0 0.0f0 0.0f0]

    G = [-sys.mt*sys.g*sin(sys.γ), 
        sys.mt*sys.g*cos(sys.γ), 
        sys.m2*sys.g*sys.l2*sin(ϕ - sys.γ),
        0.0f0]

    return B*control(z, param; expert=expert) - C*v - G
end

function impactMap(sys::RimlessWheel, ϕ)

    det = sys.I1*sys.I2 + sys.I1*sys.m2*sys.l2^2 + sys.I2*sys.mt*sys.l1^2 + 
            sys.m2*sys.l1^2*sys.l2^2*(sys.m1 + sys.m2*(sin(sys.α - ϕ))^2)

    ξ1 = 1/det * ((sys.I1*sys.I2 + sys.I1*sys.m2*sys.l2^2) + 
            (sys.I2*sys.mt*sys.l1^2 + sys.m2*sys.l1^2*sys.l2^2*(sys.m1 +0.5* sys.m2))*cos(2sys.α) -
            1/2*sys.m2^2*sys.l1^2*sys.l2^2*cos(2ϕ))

    ξ2 = 1/det *(sys.m2*sys.l1*sys.l2*(sys.I1*(cos(sys.α - ϕ) - cos(sys.α + ϕ)) + 
            sys.mt*sys.l1^2*(cos(2sys.α)*cos(sys.α - ϕ) - cos(sys.α + ϕ))) )

    return [ξ1 0;
            ξ2 1]

end

function comparePostImpact(lcp, x, param; Δt = 0.001)
    
    ϕ = x[3]
    θdot, ϕdot = [x[8], x[7]]
    postImpactMap = impactMap(sys, ϕ) * [θdot, ϕdot]

    gn, γn, γt, M, h, Wn, Wt = sysAttributes(lcp, x, param)
    λn, λt, λR  = solveLcp(lcp, gn, γn, γt, M, h, Wn, Wt; Δt=Δt)

    qM, uA = lcp.sys(x)
    uE = M\((Wn - Wt*diagm(0 => lcp.sys.μ))*λn + Wt*λR + h*Δt) + uA
    θdotlcp, ϕdotLcp  = [uE[4], uE[3]]

    println("[θdot-, ϕdot-] = ", [θdot, ϕdot])
    println("Impact map = ", postImpactMap)
    println("LCP = ", [θdotlcp, ϕdotLcp])

end

function wn(sys::RimlessWheel, z::Vector{T}) where {T<:Real}

    θ        = z[4]
    k_range  = range(1, stop=sys.k, step=1)
    Wn       = zeros(T, 4, sys.k)
    Wn[2, :] = ones(T, sys.k)
    Wn[4, :] = -sys.l1 .* sin.(θ .+ 2.0*sys.α.*k_range)

   return Wn 

end

function wt(sys::RimlessWheel, z::Vector{T}) where {T<:Real}

    θ        = z[4]
    k_range  = range(1, stop=sys.k, step=1)
    Wt       = zeros(T, 4, sys.k)
    Wt[1, :] = ones(T, sys.k)
    Wt[4, :] = -sys.l1 .* cos.(θ .+ 2.0*sys.α.*k_range)

    return Wt
end

function plots(sys, Z)
    fig1 = plt.figure()
    fig1.clf()
    subplot(2, 1, 1)
    plot(getindex.(Z, 3), getindex.(Z, 7))
    ylabel(L"\dot{\phi} [rad/s]", fontsize=15)
    subplot(2, 1, 2)
    plot(getindex.(Z, 4), getindex.(Z, 8))
    ylabel(L"\dot{\theta} [rad/s]", fontsize=15)
    fig2 = plt.figure()
    fig2.clf()
    θ = spokeInContact(sys, getindex.(Z, 4), getindex.(Z, 8))
    plot(θ, getindex.(Z, 8))
end

function plots(sys, Z, t, Λn, Λt)
    fig1 = plt.figure()
    fig1.clf()
    subplot(2, 3, 1)
    plot(t, getindex.(Z, 1))
    ylabel("x [m]", fontsize=15)
    subplot(2, 3, 2)
    plot(t, getindex.(Z, 2))
    ylabel("y [m]", fontsize=15)
    subplot(2, 3, 3)
    plot(t, getindex.(Z, 5))
    ylabel("vx [m/s]", fontsize=15)
    subplot(2, 3, 4)
    plot(t, getindex.(Z, 6))
    ylabel("vy [m/s]", fontsize=15)
    subplot(2, 3, 5)
    plot(getindex.(Z, 3), getindex.(Z, 7))
    ylabel(L"\dot{\phi} [rad/s]", fontsize=15)
    subplot(2, 3, 6)
    plot(getindex.(Z, 4), getindex.(Z, 8))
    ylabel(L"\dot{\theta} [rad/s]", fontsize=15)

    fig2 = plt.figure()
    fig2.clf()
    subplot(2, 1, 1)
    for i in 1:length(Λn[1])
        plot(t, getindex.(Λn, i))
    end
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{n1} [N]$", fontsize=15)

    subplot(2, 1, 2)
    for i in 1:length(Λt[1])
        plot(t, getindex.(Λt, i))
    end
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{t1} [N]$", fontsize=15)

    fig3 = plt.figure()
    fig3.clf()
    θ = spokeInContact(sys, getindex.(Z, 4), getindex.(Z, 8))
    plot(θ, getindex.(Z, 8))

end

function spokeInContact(sys, θ, θdot)
    θ_new = deepcopy(θ)

    for i in 2:length(θ_new)
        #if velocity jumps, wrap. Checking θ causes the angle to jump when the velocity has not
        if (abs(θdot[i] - θdot[i-1]) > 0.1) && abs(round(abs.(θ_new[i]/sys.α)) - abs.(θ_new[i]/sys.α)) < 0.03       #instead check if new contact occurs
            θ_new[i:end] .= θ_new[i:end] .+ 2*sys.α
        elseif (abs(θdot[i] - θdot[i-1]) > 0.1)  && abs(round(abs.(θ_new[i]/sys.α)) - abs.(θ_new[i]/sys.α)) < 0.03      
            θ_new[i:end] .= θ_new[i:end] .- 2sys.α
        end
    end
    return θ_new
end