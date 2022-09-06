using LaTeXStrings, PyPlot
using JuMP, PATHSolver, Mosek, MosekTools, DiffOpt
using LinearAlgebra
import ChainRulesCore, ForwardDiff
using PyCall
@pyimport matplotlib.animation as anim
using Flux, Printf


include("dynamics.jl")
include("contactMap.jl")

Δt = 0.001; totalTimeStep = 1000
θ0 = zeros(Float64, 2)

sys  = RimlessWheel(Float64, Δt, totalTimeStep)
cm  = ContactMap(Float64, sys, θ0)
x0 = [0.2, 0.1, -2.0, 0.0]

# X, t = fulltimestep(cm, x0, Float64[])
# plots(X, t)

function loss(x0::Vector{T}, S, param) where {T<:Real}
    X, t = fulltimestep(cm, x0, param)

    Q = diagm(0 => [2.0, 1.0, 0.5, 0.5])
    l = 5.0/length(X) * sum(map((x, s) -> dot(x-s, Q*(x-s)), X, S))
    return l
end

function matchExpert(cm::ContactMap, x0::Vector{T}) where {T<:Real}

    S, t = fulltimestep(cm, x0, [20.0, 5.0])

    param = [50.0, 10.0]
    l(θ) = loss(x0, S, θ)
    opt = Adam(0.1)

    while l(param) > 0.001

        grad = ForwardDiff.gradient(l, param)
        Flux.Optimise.update!(opt, param, grad)

        println("loss = ", l(param), " grad = ", grad, " θ = ", param)
    end

    return param
end

