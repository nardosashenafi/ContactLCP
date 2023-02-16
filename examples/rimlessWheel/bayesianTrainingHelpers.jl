using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra
using ForwardDiff
using Flux

include("dynamics.jl")
include("../../src/lcp.jl")
include("../../src/solver.jl")

sys            = RimlessWheel()
lcp            = Lcp(Float32, sys)

oneStep(x, γi, θ; kwargs...) = oneTimeStep(lcp, x, γi, θ; kwargs...)

function trajectory(x0, sysParam, controlParam::Vector{T}; expert=false, Δt = 0.0005f0, totalTimeStep = 1000) where {T<:Real}

    X       = Vector{Vector{T}}(undef, totalTimeStep)
    x       = deepcopy(x0)
 
    for i in 1:totalTimeStep
        x       = oneStep(x, sysParam, controlParam; Δt=Δt, expert=expert)
        X[i]    = x
    end
    # isStumbling(X) ? println("We got stumbling! Beware of gradient jumps") : nothing

    return X
end

function sampleSystemParam()
    γi = rand(0.0f0:0.01:30.0f0*pi/180f0)
    return γi
end

function sampleInitialStates(sysParam, controlParam::Vector{T}, sampleNum; α=α, totalTime=1000) where {T<:Real}

    sampleTrajectories = Vector{Vector{Vector{T}}}()

    for i in 1:2
        x0 = Float32.(initialState(rand(pi-α:0.05:pi+α), 
                                    rand(-3.0:0.05:-0.5), 
                                    0.0, 
                                    rand(-1.0:0.1:1.0)); γi = sysParam)


        S = trajectory(x0, sysParam, controlParam; totalTimeStep=totalTime)
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

    return X0

end