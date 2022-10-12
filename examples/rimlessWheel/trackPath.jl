using ContactLCP
using LaTeXStrings, PyPlot
using JuMP, LinearAlgebra
using ForwardDiff

include("dynamics.jl")
# include("../../src/lcp.jl")
# include("../../src/solver.jl")

Δt = 0.001f0; totalTimeStep = 1500
sys  = RimlessWheel(Float32)
x0 = initialState(sys, 0.0f0, -0.5f0, 0.0f0, 0.0f0)

lcp  = ContactLCP.Lcp(Float32, sys)

function loss(lcp, θ, x0, Sθd, λθd)
    Sθ, _, λθ, _ = trajectory(lcp, x0, θ)
    return 50.0f0/length(Sθ) * ( 0.5f0*dot(Sθd - Sθ , Sθd - Sθ) + 0.5f0*dot(λθd - λθ, λθd - λθ) )
end

function trackExpert(lcp::ContactLCP.Lcp, x0::Vector{T}) where {T<:Real}

    param_expert   = Float32[30.0, 5.0]
    param          = Float32[10.0, 1.0]
    η              = 0.2

    # x0 = [0.0, rand(range(0.2, 0.5, step=0.01)), 0.1, -0.1]
    Sθd, _, λθd, _ = ContactLCP.fulltimestep(lcp, x0, param_expert; Δt = 0.001f0, totalTimeStep = 1500)
    l1 = Inf
    for i in 1:20
        l(θ)  = loss(lcp, θ, x0, Sθd, λθd)
        lg    = ForwardDiff.gradient(l, param)
        θm    .-= η*lg
        l1    = l(param)

        println("loss = ", l1, " | grad = ", lg, " | param = ", param)
    end
    println("param_expert= ", param_expert, " param_learned = ", param )
    return param
end

# X, t, Λn, Λt = trajectory(lcp, x0, θ0)
# plots(X, t, Λn, Λt)
# solveM(lcp, x0)

function checkGradient(lcp::ContactLCP.Lcp, x0::Vector{T}) where {T<:Real}

    param_expert   = Float32[30.0, 5.0]
    θ1             = Float32[10.0, 1.0]
    θ2             = θ1 .+ 0.001

    Sθd, _, λθd, _ = ContactLCP.fulltimestep(lcp, x0, param_expert; Δt = 0.001f0, totalTimeStep = 1500)
    l(θ)    = loss(lcp, θ, x0, Sθd, λθd)

    l1    = l(θ1)
    l2    = l(θ2)

    lg      = ForwardDiff.gradient(l, θ1)
    # l, back = Zygote.pullback(l, θ1)
    # lg = back(1)

    fd_l = (l2 - l1) / (θ2 - θ1)
    error_grad = abs.(fd_l - lg)

    println("error_grad = ", error_grad, " | fd_l = ", fd_l, " | lg = ", lg)
end
