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
    α               ::T
    γ               ::T
    ϵn              ::Vector{T}
    ϵt              ::Vector{T}
    μ               ::Vector{T}
    x0              ::Vector{T}
    contactIndex    ::Vector{T}
    gThreshold      ::T
    Δt              ::T
    totalTimeStep   ::Int64
    stateLength     ::Int64

    function RimlessWheel(T, Δt, totalTimeStep)

        m1              = 2.32     
        m2              = 4.194 
        I1              = 0.0784160
        I2              = 0.0380256
        mt              = m1 + m2
        l1              = 0.26
        l2              = 0.05
        g               = T(9.81)
        α               = T(360.0/10.0/2.0 * pi/180.0)
        γ               = T(0.0*pi/180.0)
        ϵn              = T.(0.0*ones(1))
        ϵt              = T.(-0.5*ones(1))
        μ               = T.(0.2*ones(1))
        x0              = T.([0.0, 0.0, -0.1, 0.0])     #θ, ϕ, θdot, ϕdot
        contactIndex    = zeros(T, 1)
        gThreshold      = T.(0.001)
        Δt              = Δt
        totalTimeStep   = totalTimeStep
        stateLength     = length(x0)

        new{T}(m1, m2, I1, I2, mt, l1, l2, g, α, γ, ϵn, ϵt, μ, x0, contactIndex, gThreshold, Δt, totalTimeStep, stateLength)
    end

end

#returns the attributes need to model contact

function (sys::RimlessWheel)(x::Vector{T}, θ; limitcycle=false) where {T<:Real}
    gn  = gap(sys, x)  
    γn  = vnormal(sys, x)
    γt  = vtang(sys, x)
    M   = massMatrix(sys, x)
    h   = genForces(sys, x, θ;  limitcycle=limitcycle)
    Wn  = wn(sys, x)
    Wt  = wt(sys, x)

    return gn, γn, γt, M, h, Wn, Wt
end

#parses out the poses and velocities from the state vector
function (sys::RimlessWheel)(x::Vector{T}) where {T<:Real}
    q = x[1:2]
    u = x[3:4]
    return q, u
end

function setCoefficients(sys::RimlessWheel, ϵn, ϵt, μ)
    sys.ϵn  = ϵn
    sys.ϵt  = ϵt
    sys.μ   = μ
end

function setInitial(sys::RimlessWheel, x)
    sys.x0 = x
end

function gap(sys::RimlessWheel, x)
    return [sys.l1*cos(x[1]) - sys.l1*cos(2sys.α - abs(x[1]))] #l1cos(θ) - l1*cos(2α - |θ|)
    # return [sys.l1*(1 - cos(x[2]))]
    # return [0.0]
end

function vnormal(sys::RimlessWheel, x)
    Wn = wn(sys, x)
    _, u = sys(x)
    return Wn'*u
end

function vtang(sys::RimlessWheel, x)
    Wt = wt(sys, x)
    _, u = sys(x)
    return Wt'*u
end

function massMatrix(sys, x)
    θ, ϕ = x[1:2]
    return [sys.I1 + sys.mt*sys.l1^2 -sys.m2*sys.l1*sys.l2*cos(θ - ϕ);
            -sys.m2*sys.l1*sys.l2*cos(θ - ϕ) sys.I2 + sys.m2*sys.l2^2]
end

inputLayer(x) = [cos(x[1]), sin(x[1]), cos(x[2]), sin(x[2]), x[3], x[4]]

function control(x, θp; limitcycle=false)
   
    if limitcycle #working control for limit cycle
        return -θp[1]*(x[2]-0.325) - θp[2]*x[4]
    else
        return 0.0
        # @assert length(θp) == 6 + DiffEqFlux.paramlength(Hd)
        # y = MLBasedESC.controller(npbc, inputLayer(x), θp)
        # y = unn(inputLayer(x), θp)[1]
        # return clamp(y, -satu, satu)
    end
end

function genForces(sys, x, param; limitcycle=limitcycle)

    q, u = sys(x)
    θ, ϕ = q[1:2]
    #h = Bu - C qdot - G
    B = [-1, 1] 
    C = sys.m2*sys.l1*sys.l2*sin(θ - ϕ)*[0.0 -u[2];
                                        u[1] 0.0]

    G = sys.g*[-sys.mt*sys.l1*sin(θ - sys.γ),
                sys.m2*sys.l2*sin(ϕ - sys.γ)]

    return B*control(x, param; limitcycle=limitcycle) - C*u - G
end

function impactMap(sys, ϕ)

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

function wn(sys, x)
    sgn = sign(x[1])
    # return permutedims(hcat([0.0, -sys.l1*sin(x[1]) + sys.l1*sin(2sys.α + x[1])]...))
    return permutedims(hcat([-sys.l1*sin(x[1]) - sgn*sys.l1*sin(2sys.α - abs(x[1])), 0.0]...))
    # return permutedims(hcat([sys.l1*sin(x[1]), 0.0]...))
end

function wt(sys, x)
    return permutedims(hcat([sys.l1*cos(x[2]), 0.0]...))
end
