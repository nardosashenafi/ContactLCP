using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra
using ForwardDiff
using Flux

include("../dynamics.jl")
include("../../../src/lcp.jl")
include("../../../src/solver.jl")

sys  = RimlessWheel()
lcp  = Lcp(Float32, sys)
Δt   = 0.0003f0

function stateAndForcesWithNoise(lcp::Lcp, x, sysParam, controlParam::AbstractArray{T}; Δt = 0.001f0, kwargs...) where {T<:Real}

    qA, uA  = lcp.sys(x)
    qM      = qA + 0.5f0*Δt*uA

    x_mid   = vcat(qM, uA)
    gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold = sysAttributes(lcp, x_mid, sysParam, controlParam; kwargs...)
    λn, λt, λR  = solveLcp(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold, x_mid; Δt=Δt)

    ####complete integration
    uE = M\((Wn - Wt*diagm(0 => μ))*λn + Wt*λR + h*Δt) + uA
    qE = qM + 0.5f0*Δt*uE

    return vcat(qE,uE), λn, λt
end

function oneTimeStepWithNoise(lcp::Lcp, x, sysParam, controlParam::AbstractArray{T}; Δt = 0.001f0, kwargs...) where {T<:Real}
    x2, _, _ = stateAndForcesWithNoise(lcp, x, sysParam, controlParam; Δt = Δt, kwargs...)
    return x2
end

oneStep(x, sysParam, controlParam; kwargs...) = oneTimeStepWithNoise(lcp, x, sysParam, controlParam; kwargs...)

function trajectory(x0, r, controlParam::Vector{T}; expert=false, Δt = Δt, totalTimeStep = 1000) where {T<:Real}

    X           = Vector{Vector{T}}(undef, totalTimeStep)
    obstacles   = Vector{Vector{Vector{T}}}(undef, totalTimeStep)
    x           = deepcopy(x0)
    sysParam    = createUnevenTerrain(x, r)
    θi          = x[4]

    for i in 1:totalTimeStep

        if abs(x[4] - θi) > 2pi - 2α
            sysParam = createUnevenTerrain(x, r)
            θi = x[4]
        end
        x       = oneStep(x, sysParam, controlParam; Δt=Δt, expert=expert)
        X[i]    = x
        obstacles[i] = sysParam
    end

    return X, obstacles
end

function sampleSystemParam()
    γi = rand(0.0f0:0.01:30.0f0*pi/180f0)
    return γi
end

function isStumbling(x)
    x[8] >= -0.01
end

function sampleInitialStates(controlParam::Vector{T}, sampleNum; α=α, totalTime=1000) where {T<:Real}

    sampleTrajectories = Vector{Vector{Vector{T}}}()
    rmax = 0.06f0

    w     = rand(getq(controlParam))
    x0, r = initialStateWithBumps(rand(pi-α:0.05f0:pi+α), 
                                rand(-3.0f0:0.05f0:-0.5f0), 
                                0.0f0, 
                                rand(-1.0f0:0.1f0:1.0f0), rmax)


    S, _ = trajectory(x0, r, w; totalTimeStep=totalTime)
    push!(sampleTrajectories, S)

    #sample 10 initial states from the long trajectories
    X0 = Vector{Vector{T}}()
    R = Vector{Vector{T}}()

    for i in 1:sampleNum
        if rand() < 0.5 
            push!(X0, rand(rand(sampleTrajectories)))
            push!(R, r)
        else
            x0, r1 = initialStateWithBumps(rand(pi-α:0.05f0:pi+α), 
                                        rand(-3.0f0:0.05f0:-0.5f0), 
                                        0.0f0, 
                                        rand(-1.0f0:0.1f0:1.0f0), rmax)
            push!(X0, x0)
            push!(R, r1)
        end 
    end
    #extract stumbling
    for i in eachindex(X0)
        if isStumbling(X0[i])
            X0[i], R[i] = Float32.(initialState(rand(pi-α:0.05f0:pi+α), 
                        rand(-3.0f0:0.05f0:-0.5f0), 
                        rand(-pi/4.0f0:0.05f0:pi/4.0f0), 
                        rand(-1.0f0:0.1f0:1.0f0)), rmax)
        end
    end
    return X0, R

end

function plots(Z, fig1)
    PyPlot.figure(1)
    fig1.clf()
    subplot(2, 2, 1)
    plot(getindex.(Z, 3), getindex.(Z, 7))
    scatter(Z[end][3], Z[end][7])
    ylabel(L"\dot{\phi} [rad/s]", fontsize=15)
    subplot(2, 2, 2)
    # θ, impactIndex = spokeInContact(Z)
    plot(getindex.(Z, 4), getindex.(Z, 8))
    scatter(Z[end][4], Z[end][8])
    ylabel(L"\dot{\theta} [rad/s]", fontsize=15)
    subplot(2, 2, 3)
    plot(getindex.(Z, 5))
    ylabel("vx [m/s]", fontsize=15)
    println("Average hip speed = ", mean(getindex.(Z, 5)))
end