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

function loss(lcp, param, x0; totalTime = 500)
    S, λ, _ = rwTrajectory(lcp, x0, param; totalTimeStep=totalTime)
    θ, θdot, ϕ, ϕdot = getindex.(S, 4), getindex.(S, 8), getindex.(S, 3), getindex.(S, 7)
    Sd, λd, _ = rwTrajectory(lcp, x0, param_expert; totalTimeStep=totalTime); 
    θd, θdotd, ϕd, ϕdotd = getindex.(Sd, 4), getindex.(Sd, 8), getindex.(Sd, 3), getindex.(Sd, 7)

    return 500.0f0/length(S) * ( 
            2.0f0*dot(θd - θ , θd - θ) +
            2.0f0*dot(θdotd - θdot , θdotd - θdot) + 
            2.0f0*dot(ϕd - ϕ , ϕd - ϕ) + 
            2.0f0*dot(ϕdotd - ϕdot , ϕdotd - ϕdot) +  
            0.0f0*dot(λd - λ, λd - λ) )
end

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

    return X, Λn, Λt
end

function isStumbling(x)
    # length(x) <= 50
    any(getindex.(x, 8) .> 0.0) 
    # false
end

function sampleInitialStates(lcp::Lcp, param::Vector{T}; totalTime=1000) where {T<:Real}

    sampleTrajectories = Vector{Vector{Vector{T}}}()

    #generate 5 long trajectories
    while length(sampleTrajectories) < 2

        x0 = initialState(sys, 0.0f0, -0.5f0, 0.0f0, 0.0f0)
        S, λ, _ = rwTrajectory(lcp, x0, param; totalTimeStep=totalTime)
        isStumbling(S) ? println("We got stumbling! Beware of gradient jumps") : nothing
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

function trackExpert(lcp::Lcp, x0::Vector{T}) where {T<:Real}
    param       = Float32[25.0, 2.0]
    opt         = Adam(0.01)

    l1          = Inf
    counter     = 0
    X0          = Vector{Vector{T}}()

    for i in 1:400
        while isempty(X0)
            X0 = sampleInitialStates(lcp, param; totalTime=1000)
        end
        for xi in X0
            l(θ)  = loss(lcp, θ, xi; totalTime = 500)
            lg    = ForwardDiff.gradient(l, param)

            l1    = l(param)
            if counter > 20
                println("loss = ", l1, " | grad = ", lg, " | param = ", param)
                counter = 0
            end
            counter += 1

            Flux.update!(opt, param, lg)
        end
    end
    println("param_expert= ", param_expert, " param_learned = ", param )
    return param
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


