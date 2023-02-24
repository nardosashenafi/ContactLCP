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
    uE = M\((Wn_noise - Wt*diagm(0 => μ))*λn + Wt_noise*λR + h*Δt) + uA
    qE = qM + 0.5f0*Δt*uE

    return vcat(qE,uE), λn, λt
end

function oneTimeStepWithNoise(lcp::Lcp, x, sysParam, controlParam::AbstractArray{T}; Δt = 0.001f0, kwargs...) where {T<:Real}
    x2, _, _ = stateAndForcesWithNoise(lcp, x, sysParam, controlParam; Δt = Δt, kwargs...)
    return x2
end

oneStep(x, sysParam, controlParam; kwargs...) = oneTimeStepWithNoise(lcp, x, sysParam, controlParam; kwargs...)

function bumpHeight(x)

    x_floor, y_floor = [0.0, 0.0]
    gn, _, _ = gap(x)
    hs = hangingSpoke(x, gn)
    xs = x[1] - l1*sin(θ + 2*α*hs)
    gi, _, _ = gapHangingSpoke(x, hs, x_floor, y_floor) 

end

function trajectory(x0, controlParam::Vector{T}; expert=false, Δt = Δt, totalTimeStep = 1000) where {T<:Real}

    X       = Vector{Vector{T}}(undef, totalTimeStep)
    x       = deepcopy(x0)
    sysParam = bumpHeight(x)

    for i in 1:totalTimeStep
        x       = oneStep(x, sysParam, controlParam; Δt=Δt, expert=expert)
        X[i]    = x
    end

    return X
end

function sampleSystemParam()
    γi = rand(0.0f0:0.01:30.0f0*pi/180f0)
    return γi
end

function isStumbling(x)
    x[8] >= -0.01
end

function sampleInitialStates(sysParam, controlParam::Vector{T}, sampleNum; α=α, totalTime=1000) where {T<:Real}

    sampleTrajectories = Vector{Vector{Vector{T}}}()

    for i in 1:2
        w  = rand(getq(controlParam))
        x0 = Float32.(initialState(rand(pi-α:0.05:pi+α), 
                                    rand(-3.0:0.05:-0.5), 
                                    0.0, 
                                    rand(-1.0:0.1:1.0)))


        S = trajectory(x0, w; totalTimeStep=totalTime)
        push!(sampleTrajectories, S)
    end

    #sample 10 initial states from the long trajectories
    X0 = Vector{Vector{T}}()

    for i in 1:sampleNum
        if rand() < 0.5 
            push!(X0, rand(rand(sampleTrajectories)))
        else
            x0 = Float32.(initialState(rand(pi-α:0.05:pi+α), 
                                        rand(-3.0:0.05:-0.5), 
                                        0.0, 
                                        rand(-1.0:0.1:1.0)) )
            push!(X0, x0)
        end 
    end
    #extract stumbling
    for i in eachindex(X0)
        if isStumbling(X0[i])
            X0[i] = Float32.(initialState(rand(pi-α:0.05:pi+α), 
                        rand(-3.0:0.05:-0.5), 
                        rand(-pi/4:0.05:pi/4), 
                        rand(-1.0:0.1:1.0)) )
        end
    end
    return X0

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