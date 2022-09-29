export BouncingBall

mutable struct BouncingBall{T} 
    m               ::T 
    r               ::T
    g               ::T 
    ϵn              ::Vector{T}
    ϵt              ::Vector{T}
    μ               ::Vector{T}
    x0              ::Vector{T}
    contactIndex    ::Vector{T}
    gThreshold      ::T
    Δt              ::T
    totalTimeStep   ::Int64
    stateLength     ::Int64

    function BouncingBall(T, Δt, totalTimeStep)

        m               = T(0.21)
        r               = T(0.1)
        g               = T(9.81)
        ϵn              = T.(0.5*ones(1))
        ϵt              = T.(-0.5*ones(1))
        μ               = T.(0.2*ones(1))
        x0              = T.([0.0, 0.5, 0.1, -0.1])     #xpos, ypos, xdot, ydot
        contactIndex    = zeros(T, 1)
        gThreshold      = T.(0.001)
        Δt              = Δt
        totalTimeStep   = totalTimeStep
        stateLength     = length(x0)

        new{T}(m, r, g, ϵn, ϵt, μ, x0, contactIndex, gThreshold, Δt, totalTimeStep, stateLength)
    end

end

#returns the attributes need to model contact
function (sys::BouncingBall)(x::Vector{T}, θ::Vector{T}) where {T<:Real}
    gn  = gap(sys, x)
    γn  = vnormal(sys, x)
    γt  = vtang(sys, x)
    M   = massMatrix(sys, x, θ)
    h   = genForces(sys, x, θ)
    Wn  = wn(sys, x)
    Wt  = wt(sys, x)

    return gn, γn, γt, M, h, Wn, Wt
end

#parses out the poses and velocities from the state vector
function (sys::BouncingBall)(x::Vector{T}) where {T<:Real}
    q = x[1:2]
    u = x[3:4]
    return q, u
end

# returns the perturbation or derivatives of the objective and constraints of the LCP optimization with respect to 
# the learned parameters. This is passed to the ContactLCP package through the function getPerturbations()
function (sys::BouncingBall)(model, x::Vector{T}, θ::Vector{T}, A::Matrix{T}, b::Vector{T}, sol) where {T<:Real}

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

function setCoefficients(sys::BouncingBall, ϵn, ϵt, μ)
    sys.ϵn = ϵn
    sys.ϵt = ϵt
    sys.μ= μ
end

function setInitial(sys::BouncingBall, x)
    sys.x0 = x
end

function setSysParams(sys::BouncingBall, params)
    sys.m = params[1]
    sys.r = params[2]
end

function gap(sys::BouncingBall, x)
    return [x[2] - sys.r]
end

function vnormal(sys::BouncingBall, x)
    Wn = wn(sys, x)
    _, u = sys(x)
    return Wn'*u
end

function vtang(sys::BouncingBall, x)
    Wt = wt(sys, x)
    _, u = sys(x)
    return Wt'*u
end

function massMatrix(sys, x, θm)
    return [θm[1] 0.0;
            0.0 θm[1]]
end

function genForces(sys, x, θm)
    return [0.0, -θm[1]*sys.g]
end

function wn(sys, x)
    return permutedims(hcat([0.0, 1.0]...))
end

function wt(sys, x)
    return permutedims(hcat([1.0, 0.0]...))
end

function plots(Z, t, Λn)
    fig1 = plt.figure()
    fig1.clf()
    subplot(2, 2, 1)
    plot(t, getindex.(Z, 1))
    ylabel("x [m]", fontsize=15)
    subplot(2, 2, 2)
    plot(t, getindex.(Z, 2))
    ylabel("y [m]", fontsize=15)
    subplot(2, 2, 3)
    plot(t, getindex.(Z, 3))
    ylabel("vx [m/s]", fontsize=15)
    subplot(2, 2, 4)
    plot(t, getindex.(Z, 4))
    ylabel("vy [m/s]", fontsize=15)

    fig2 = plt.figure()
    fig2.clf()
    subplot(2, 1, 1)
    plot(t[2:end], getindex.(Λn, 1))
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{n1} [N]$", fontsize=15)

end
