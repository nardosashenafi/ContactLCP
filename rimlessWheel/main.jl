using LaTeXStrings, PyPlot
using JuMP, PATHSolver, Mosek, MosekTools, DiffOpt
using LinearAlgebra
import ChainRulesCore, ForwardDiff
using PyCall
@pyimport matplotlib.animation as anim
using Flux, Printf
using Interpolations
using Lux
using Statistics

include("dynamics.jl")
include("contactMap.jl")

Δt = 0.001; totalTimeStep = 500
# unn     = Lux.Chain(Lux.Dense(4, 5, relu), Lux.Dense(5,1))
# controlparam = rand(Lux.parameterlength(unn))
controlparam = [5.0, 1.0, 0.2]
θ0  = zeros(Float64, length(controlparam))

sys = RimlessWheel(Float64, Δt, totalTimeStep)
cm  = ContactMap(Float64, sys, θ0)
x0  = [0.2, cm.sys.γ, -2.0, 0.0]
replayBuffer = []

# X, t = fulltimestep(cm, x0, Float64[])
# plots(X, t)

function interpolate_τ(limitcycle_x, limitcycle_t)
    limitcycle_θ    =  getindex.(limitcycle_x, 1)
    τ               = Interpolations.extrapolate(Interpolations.interpolate((-limitcycle_θ,), limitcycle_t, Gridded(Linear())), Line())
    return τ
end

function interpolate_θdot(limitcycle_x, limitcycle_t)
    limitcycle_θdot =  getindex.(limitcycle_x, 3)
    θdot_star       = Interpolations.extrapolate(Interpolations.interpolate((limitcycle_t,), limitcycle_θdot, Gridded(Linear())), Line())
    return θdot_star
end

function isStumbling(x)
    # length(x) <= 50
    (any(getindex.(x, 3) .> 0.0) || length(x) <= 2)     
    # false
end

function perpLoss(x0, τ, θdotstar, param::Vector{T}) where {T<:Real}

    X, t        = fulltimestep(cm, x0, param; timeSteps=500)
    loss        = 0.0
    Z           = Vector{Vector{Vector{T}}}()

    for x in X
        if !isStumbling(x)
            push!(Z, x)
        end
    end

    for x in Z
        θ           = getindex.(x, 1)
        θdot        = getindex.(x, 3)
        ϕdot        = getindex.(x, 4)
        τi          = τ.(-θ)
        θdotstar_i  = θdotstar.(τi)
        loss        += 1.0/length(θ)*(sum(map((xs, s) -> dot(xs-s, xs-s), θdot, θdotstar_i)) +
                        0.0*dot(ϕdot, ϕdot))

    end

    return loss/length(Z)
end


function hipSpeedLoss(cm, x0, param::Vector{T}) where {T<:Real}
    xd_dot = 0.5
    X, t   = fulltimestep(cm, x0, param; timeSteps=4000)
    Z      = Vector{Vector{Vector{T}}}()
    for x in X
        if !isStumbling(x)
            push!(Z, x)
        end
    end

    loss = 0.0
    for x in Z
        θ       = getindex.(x, 1)
        θdot    = getindex.(x, 3)
        ϕdot    = getindex.(x, 4)
        l1      = xd_dot .+ cm.sys.l1 .* cos.(θ) .* θdot 
        loss    += 1/length(θ)*dot(l1, l1) 
    end

    return 100.0/length(Z)*loss
end

function sampleInitialStates(x0, param::Vector{T}) where {T<:Real}
    X, t = fulltimestep(cm, x0, param; timeSteps=5000)
    Z    = Vector{Vector{Vector{T}}}()

    for x in X
        if !isStumbling(x)
            push!(Z, x)
        end
    end

    X0 = Vector{Vector{T}}()
    sampleNum = 50
    for i in 1:sampleNum
        push!(X0, rand(rand(Z)))
    end
    return X0

end

function controlToHipSpeed(cm::ContactMap)

    counter     = 0
    param       = deepcopy(controlparam)
    opt         = Adam(0.005)
    fig1        = plt.figure()
    loss        = Inf
    X0          = Vector{Vector{Float64}}()

    for i in 1:2000
        while isempty(X0)
            # x0 = [rand(-cm.sys.α:0.05:cm.sys.α), rand(0:0.1:pi/2), rand(-5.0:0.1:0.0), 0.0]
            X0 = sampleInitialStates(x0, param)
        end

        for xi in X0
            l1(θ)   = hipSpeedLoss(cm, xi, θ)
            grad    = ForwardDiff.gradient(l1, param)
            
            if counter > 30
                X, t    = fulltimestep(cm, xi, param; timeSteps=5000)
                Z       = reduce(vcat, X)
                fig1.clf()
                subplot(2, 1, 1)
                plot(getindex.(Z, 1), getindex.(Z, 3))
                subplot(2, 1, 2)
                plot(getindex.(Z, 2))

                # @info "loss" l1(param) "g" grad "hip speed" mean(cm.sys.l1 .* cos.(getindex.(Z, 1)) .* getindex.(Z, 3)) 
                println("loss = ", round(l1(param), digits=4),  " p = ", round.(param, digits=4), " hip speed = ", round.(mean(-cm.sys.l1 .* cos.(getindex.(Z, 1)) .* getindex.(Z, 3)), digits=4) )
                counter = 0
            end

            Flux.Optimise.update!(opt, param, grad)
            counter += 1
            # loss = l1(param)
        end
    end

end

function controlToLimitCycle(cm::ContactMap) 

    S, t        = limitCycle(cm)
    τ           = interpolate_τ(S, t)
    θdotstar    = interpolate_θdot(S, t)
    counter     = 0
    param       = deepcopy(controlparam)
    opt         = Adam(0.005)
    fig1        = plt.figure()
    loss        = Inf
    X0          = Vector{Vector{T}}()

    for i in 1:2000
        while isempty(X0)
            # x0 = [rand(-cm.sys.α:0.05:cm.sys.α), rand(0:0.1:pi/2), rand(-5.0:0.1:0.0), 0.0]
            X0 = sampleInitialStates(x0, param)
        end

        for xi in X0
            l1(θ)   = perpLoss(xi, τ, θdotstar , θ)
            grad    = ForwardDiff.gradient(l1, param)
            if counter > 10
                println("loss = ", l1(param), " θ = ", param)
                println("grad = ", grad)
                X, t    = fulltimestep(cm, xi, param; timeSteps=5000)
                Z       = reduce(vcat, X)
                fig1.clf()
                subplot(2, 1, 1)
                plot(getindex.(Z, 1), getindex.(Z, 3))
                plot(getindex.(S, 1), getindex.(S, 3))
                subplot(2, 1, 2)
                plot(getindex.(Z, 2))
                counter = 0
            end

            Flux.Optimise.update!(opt, param, grad)
            counter += 1
            loss = l1(param)
        end
    end

end


function fd(x0, τ, θdotstar, Δ, param1)
    param2 = param1 + [Δ, 0.0]
    l1 = perpLoss(x0, τ, θdotstar, param1) 
    l2 = perpLoss(x0, τ, θdotstar, param2)

    param3 = param1 + [0.0, Δ]
    l3 = perpLoss(x0, τ, θdotstar, param3)
    return [(l2 - l1) ./ (param2[1] - param1[1]), (l3 - l1) ./ (param3[2] - param1[2])]
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
