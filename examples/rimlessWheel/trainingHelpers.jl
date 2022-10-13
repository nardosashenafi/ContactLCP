# using ContactLCP
using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra
using ForwardDiff
using Flux

include("dynamics.jl")
include("../../src/lcp.jl")
include("../../src/solver.jl")

Δt = 0.001f0; totalTimeStep = 1500
sys  = RimlessWheel(Float32)
x0 = initialState(sys, 0.0f0, -0.5f0, 0.0f0, 0.0f0)
param_expert   = Float32[30.0, 5.0]
lcp  = Lcp(Float32, sys)


function rwTrajectory(lcp::Lcp, x0, param::Vector{T}; Δt = 0.001f0, totalTimeStep = 1000) where {T<:Real}

    X       = Vector{Vector{T}}(undef, totalTimeStep)
    Λn      = Vector{Vector{T}}(undef, totalTimeStep)
    Λt      = Vector{Vector{T}}(undef, totalTimeStep)
    x       = deepcopy(x0)

    for i in 1:totalTimeStep
        x, λn, λt  = oneTimeStep(lcp, x, param; Δt=Δt)
        X[i]    = x
        Λn[i]   = λn
        Λt[i]   = λt
    end
    # isStumbling(X) ? println("We got stumbling! Beware of gradient jumps") : nothing

    return X, Λn, Λt
end

function isStumbling(x)
    any(getindex.(x, 8) .> 0.0) 
end

function sampleInitialStates(lcp::Lcp, param::Vector{T}; totalTime=1000) where {T<:Real}

    sampleTrajectories = Vector{Vector{Vector{T}}}()

    #generate 5 long trajectories
    while length(sampleTrajectories) < 2

        x0 = initialState(sys, 0.0f0, -0.5f0, 0.0f0, 0.0f0)
        # x0 = Float32.(initialState(sys, rand(-lcp.sys.α+0.1:0.05:0.0), 
        #                               rand(-2.0:0.05:-0.5), 0.0f0, 0.0f0))

        S, λ, _ = rwTrajectory(lcp, x0, param; totalTimeStep=totalTime)
        # isStumbling(S) ? println("We got stumbling! Beware of gradient jumps") : nothing
        # Z, tz = extractStumbling(X, tx) 
        push!(sampleTrajectories, S)

    end

    #sample 10 initial states from the long trajectories
    X0 = Vector{Vector{T}}()
    sampleNum = 4

    for i in 1:sampleNum
        push!(X0, rand(rand(sampleTrajectories)))
    end

    return X0

end

# X, t, Λn, Λt = trajectory(lcp, x0, θ0)
# plots(X, t, Λn, Λt)
# solveM(lcp, x0)

function checkRWGradient(lcp::Lcp, x0::Vector{T}) where {T<:Real}
 

    θ1     = Float32[10.0, 2.0]
    θ2     = θ1 .+ Float32[0.0001, 0.0]

    l(θ)   = loss(lcp, θ, x0)

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


