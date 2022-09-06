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
function (sys::RimlessWheel)(x::Vector{T}, θ::Vector{T}) where {T<:Real}
    gn  = gap(sys, x)
    γn  = vnormal(sys, x)
    γt  = vtang(sys, x)
    M   = massMatrix(sys, x)
    h   = genForces(sys, x, θ)
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

# returns the perturbation or derivatives of the objective and constraints of the LCP optimization with respect to 
# the learned parameters. This is passed to the ContactLCP package through the function getPerturbations()
function (sys::RimlessWheel)(model, x::Vector{T}, θ::Vector{T}, A::Matrix{T}, b::Vector{T}, sol) where {T<:Real}

    q_new, v_new, λ_new = sol
    uA = x[3:4]
    Δcons1 = -1.0/θ[1]^2.0*index(model[:λ][1]) + 0.0
    Δcons2 = [  1.0*(index(model[:v][1]) - uA[1])
                1.0*(index(model[:v][2])) - 0.0*index(model[:λ][1]) - uA[2] + sys.g*sys.Δt     #∂cons2[2]/∂θm
             ]

    Δcons3 = nothing
    Δobj   = nothing
    if !isempty(A)
        Δobj = -1.0/θ[1]^2.0 * λ_new[1]' * A[1,1] * (ones(Float64) *index(model[:λ][1])*index(model[:λ][1])) + 0.0
    end

    return Δcons1, Δcons2, Δcons3, Δobj
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

function control(x, θ)
    return 0
end

function genForces(sys, x, param)
    q, u = sys(x)
    θ, ϕ = q[1:2]
    #h = Bu - C qdot - G
    B = [-1, 1] 
    C = sys.m2*sys.l1*sys.l2*sin(θ - ϕ)*[0.0 -u[2];
                                        u[1] 0.0]

    G = sys.g*[-sys.mt*sys.l1*sin(θ - sys.γ),
                sys.m2*sys.l2*sin(ϕ - sys.γ)]

    return B*control(x, param) - C*u - G
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

function plotRimless(sys, Z)

end

function plots(Z, t, Λn)
    fig1 = plt.figure()
    fig1.clf()
    subplot(2, 2, 1)
    plot(t, getindex.(Z, 1))
    ylabel(L"$\theta$", fontsize=15)
    subplot(2, 2, 2)
    plot(t, getindex.(Z, 2))
    ylabel(L"$\phi$", fontsize=15)
    subplot(2, 2, 3)
    plot(t, getindex.(Z, 3))
    ylabel(L"$\dot{\theta}$", fontsize=15)
    subplot(2, 2, 4)
    plot(t, getindex.(Z, 4))
    ylabel(L"$\dot{\phi}$", fontsize=15)



    fig2 = plt.figure()
    fig2.clf()
    subplot(2, 1, 1)
    plot(t[2:end], getindex.(Λn, 1))
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{n1} [N]$", fontsize=15)
    subplot(2, 1, 2)
    plot(getindex.(Z, 1), getindex.(Z, 3))
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\dot{\theta}$", fontsize=15)
    xlabel(L"$\theta$", fontsize=15)

end
