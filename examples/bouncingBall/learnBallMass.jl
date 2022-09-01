
using ContactLCP
using LaTeXStrings
using PyPlot
using JuMP

include("dynamics.jl")

Δt = 0.001; totalTimeStep = 1500
θ0 = Float64[0.2]
x0 = [0.0,0.5,0.2,-0.2]

sys  = BouncingBall(Float64, Δt, totalTimeStep)
lcp  = ContactLCP.Lcp(Float64, sys, θ0)

function loss(Sθd, Sθ, λθd, λθ)
    return 0.5*sum(abs2, λθ - λθd) 
end

function ∂l∂optsol(Sθd, Sθ, λθd, λθ)
    return (λθ - λθd)'
end

function solveM(lcp::ContactLCP.Lcp, x0::Vector{T}, stateNum) where {T<:Real}

    θm_actual   = [0.2]
    θm          = [0.9]
    η           = 0.05

    # x0 = [0.0, rand(range(0.2, 0.5, step=0.01)), 0.1, -0.1]
    S, t, Λ     = ContactLCP.fulltimestep(lcp, x0, θm_actual)
    for i in 1:10
        l, lg   = ContactLCP.lossGrad(lcp, S, Λ, θm, stateNum, loss, ∂l∂optsol)
        println("loss = ", l, " | grad = ", lg, " | θm = ", θm)
        θm      .-= η*lg
    end
    println("θ_actual= ", θm_actual, " θ_learned = ", θm )
    return θm
end

# X, t, Λn = ContactLCP.fulltimestep(lcp, x0, θ0)
# plots(X, t, Λn)
solveM(lcp, x0, 2)