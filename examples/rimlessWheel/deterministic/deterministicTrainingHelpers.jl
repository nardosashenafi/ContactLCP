# using ContactLCP
using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra
using ForwardDiff
using Flux

include("../dynamics.jl")
include("../../../src/lcp.jl")
include("../../../src/solver.jl")

sys            = RimlessWheel()
lcp            = Lcp(Float32, sys)
Δt = 0.0003f0; totalTimeStep = 1500

oneStep(x, θ; kwargs...) = oneTimeStep(lcp, x, θ; kwargs...)

function trajectory(x0, param::Vector{T}; expert=false, Δt = Δt, totalTimeStep = totalTimeStep) where {T<:Real}

    X       = Vector{Vector{T}}(undef, totalTimeStep)
    x       = deepcopy(x0)
 
    for i in 1:totalTimeStep
        x       = oneStep(x, param; Δt=Δt, expert=expert)
        X[i]    = x
    end
    # isStumbling(X) ? println("We got stumbling! Beware of gradient jumps") : nothing

    return X
end

function isStumbling(x)
    x[8] >= -0.01
end

function sampleInitialStates(param::Vector{T}, sampleNum; α=α, totalTime=1000) where {T<:Real}

    sampleTrajectories = Vector{Vector{Vector{T}}}()

    # x0 = initialState(pi-0.0f0, -0.5f0, 0.0f0, 0.0f0)
    x0 = Float32.(initialStateRandomAngle(rand(pi-α:0.05:pi+α), 
                                rand(-3.0:0.05:0.0),   #large thetadot can cause leaping
                                0.0, 
                                rand(-1.0:0.1:1.0)) )

    S = trajectory(x0, param; totalTimeStep=totalTime)
    push!(sampleTrajectories, S)

    #sample 10 initial states from the long trajectories
    X0 = Vector{Vector{T}}()

    for i in 1:sampleNum
        if rand() < 0.4 
            push!(X0, rand(rand(sampleTrajectories)))
        else
            x0 = Float32.(initialStateRandomAngle(rand(pi-α:0.05:pi+α), 
                                        rand(-3.0:0.05:0.0), 
                                        0.0, 
                                        rand(-1.0:0.1:1.0)) )
            push!(X0, x0)
        end 
    end
    #extract stumbling
    for i in eachindex(X0)
        if isStumbling(X0[i])
            X0[i] = Float32.(initialStateRandomAngle(rand(pi-α:0.05:pi+α), 
                        rand(-3.0:0.05:0.0), 
                        rand(-pi/4:0.05:pi/4), 
                        rand(-1.0:0.1:1.0)) )
        end
    end
    return X0

end

# X, t, Λn, Λt = trajectory(lcp, x0, θ0)
# plots(X, t, Λn, Λt)
# solveM(lcp, x0)

function checkRWGradient(x0::Vector{T}, ps) where {T<:Real}
 

    θ1     = ps
    θ2     = deepcopy(θ1)
    θ2[1]  = θ2[1] + 0.0001

    l(θ)   = trajLoss([x0], θ; totalTime=500)

    l1    = l(θ1)
    l2    = l(θ2)
    fd_l  = (l2 - l1) / (θ2[1] - θ1[1])

    lg      = ForwardDiff.gradient(l, θ1)
    # l, back = Zygote.pullback(l, θ1)
    # lg = back(1)

    perc_error = 100.0 * abs.(fd_l - lg[1])/fd_l
    #stumbling trajectories make finite diff sensitive 
    println("perc_error = ", perc_error, " | fd_l = ", fd_l, " | lg = ", lg[1])
end


