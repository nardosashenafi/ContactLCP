using LaTeXStrings, PyPlot
# using JuMP, PATHSolver, Mosek, MosekTools, DiffOpt
using LinearAlgebra
# import ChainRulesCore
using ForwardDiff, DiffEqFlux
using PyCall
@pyimport matplotlib.animation as anim
using Flux, Printf
using Interpolations
# using Lux, Random, Optimisers, Zygote, AbstractDifferentiation
using Statistics


Δt              = 0.001 
totalTimeStep   = 500
# unn             = Lux.Chain(x -> [cos(x[1]), sin(x[1]), cos(x[2]), sin(x[2]), x[3], x[4]], 
#                     Lux.Dense(6, 3, elu), Lux.Dense(3,1))
customLayer(x, p) = [cos(x[1]), sin(x[1]), cos(x[2]), sin(x[2]), x[3], x[4]]
jacobian(::typeof(customLayer), x) = [-sin(x[1]) 0.0 0.0 0.0;
                                       cos(x[1]) 0.0 0.0 0.0;
                                       0.0 -sin(x[2]) 0.0 0.0;
                                       0.0 cos(x[2]) 0.0 0.0;
                                       0.0 0.0 1.0 0.0;
                                       0.0 0.0 0.0 1.0]
                     
# _applychain(fs::InputFastDense, x, p) = _applychain(tail(fs), first(fs)(x), p)
unn             = FastChain(customLayer, FastDense(6, 8, elu), FastDense(8, 1))
# rng             = Random.default_rng()
# ps, st          = Lux.setup(rng, unn)
# controlparam    = rand(Lux.parameterlength(unn))
ps          = 0.05*rand(DiffEqFlux.paramlength(unn))
satu        = 1.0
# controlparam = [5.0, 1.0, 0.2]

include("dynamics.jl")
include("contactMap.jl")


sys = RimlessWheel(Float64, Δt, totalTimeStep)
cm  = ContactMap(Float64, sys, ps)
x0  = [0.2, cm.sys.γ, -2.0, 0.0]

# X, t = fulltimestep(cm, x0, Float64[])
# plots(X, t)

include("learnControl.jl")