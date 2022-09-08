using LaTeXStrings, PyPlot
using JuMP, PATHSolver, Mosek, MosekTools, DiffOpt
using LinearAlgebra
import ChainRulesCore, ForwardDiff
using PyCall
@pyimport matplotlib.animation as anim
using Flux, Printf
using Interpolations

include("dynamics.jl")
include("contactMap.jl")

Δt = 0.001; totalTimeStep = 5000
θ0 = zeros(Float64, 2)

sys  = RimlessWheel(Float64, Δt, totalTimeStep)
cm  = ContactMap(Float64, sys, θ0)
x0 = [0.2, 0.1, -2.0, 0.0]

# X, t = fulltimestep(cm, x0, Float64[])
# plots(X, t)

function interpolate_τ(limitcycle_x, limitcycle_t)
    limitcycle_θ    =  getindex.(limitcycle_x, 1)
    τ               = Interpolations.extrapolate(Interpolations.interpolate((-limitcycle_θ,), limitcycle_t, Gridded(Linear())), Line())
    return τ
end

function interpolate_θdot(limitcycle_x, limitcycle_t)
    limitcycle_θdot =  getindex.(limitcycle_x, 3)
    θdot_star   = Interpolations.extrapolate(Interpolations.interpolate((limitcycle_t,), limitcycle_θdot, Gridded(Linear())), Line())
    return θdot_star
end

function perpLoss(x0, τ, θdotstar, param::Vector{T}) where {T<:Real}
    X, t        = fulltimestep(cm, x0, param)
    loss        = 0.0

    for x in X
        θ           = getindex.(x, 1)
        θdot        = getindex.(x, 3)
        ϕdot        = getindex.(x, 2)
        τi          = τ.(-θ)
        θdotstar_i  = θdotstar.(τi)
        loss        += 1.0/length(θ)*(sum(map((x, s) -> dot(x-s, x-s), θdot, θdotstar_i)) + 0.0*dot(ϕdot, ϕdot))
    end
    return loss
end

function controlToLimitCycle(cm::ContactMap) 

    S, t        = limitCycle(cm)
    τ           = interpolate_τ(S, t)
    θdotstar    = interpolate_θdot(S, t)
    counter     = 0
    param       = [5.0, 2.0]
    opt         = Adam(0.01)
    fig1 = plt.figure()
    fig1.clf()
    loss = Inf
   
    while loss > 0.001
        l1(θ)  = perpLoss(x0, τ, θdotstar , θ)
        grad = ForwardDiff.gradient(l1, param)
        Flux.Optimise.update!(opt, param, grad)
        if counter > 50
            println("loss = ", l1(param), " θ = ", param)
            println("grad = ", grad)
            X, t    = fulltimestep(cm, x0, param)
            Z       = reduce(vcat, X)
            t1      = reduce(vcat, t)
            fig1.clf()
            subplot(2, 1, 1)
            plot(getindex.(Z, 1), getindex.(Z, 3))
            plot(getindex.(S, 1), getindex.(S, 3))
            subplot(2, 1, 2)
            plot(t1, getindex.(Z, 2))
            counter = 0
        end
        counter += 1
        loss = l1(param)
    end

end

function matchExpert(cm::ContactMap, x0::Vector{T}) where {T<:Real}

    function matchloss(x0::Vector{T}, S, param) where {T<:Real}
        X, t = fulltimestep(cm, x0, param)
    
        Q = diagm(0 => [2.0, 1.0, 0.5, 0.5])
        l = 5.0/length(X) * sum(map((x, s) -> dot(x-s, Q*(x-s)), X, S))
        return l
    end

    S, t    = fulltimestep(cm, x0, [20.0, 5.0])

    param   = [50.0, 10.0]
    l1(θ)   = matchloss(x0, S, θ)
    opt     = Adam(0.1)
    counter = 0

    while l1(param) > 0.001
        grad = ForwardDiff.gradient(l1, param)
        Flux.Optimise.update!(opt, param, grad)

        if counter > 50
            println("loss = ", l1(param), " grad = ", grad, " θ = ", param)
            counter = 0
        end
        counter += 1
    end

    return param
end

