using LaTeXStrings, PyPlot
# using JuMP, PATHSolver, Mosek, MosekTools, DiffOpt
using LinearAlgebra
# import ChainRulesCore
using ForwardDiff, DiffEqFlux
using PyCall
@pyimport matplotlib.animation as anim
using Flux, Printf
using Interpolations
using Statistics
using MLBasedESC

Δt              = 0.001 
totalTimeStep   = 500

unn             = FastChain(FastDense(6, 8, elu), FastDense(8, 1))
Hd              = FastChain(FastDense(6, 12, elu), FastDense(12, 1))
N               = 6
npbc            = MLBasedESC.NeuralPBC(N, Hd)

# ps              = 0.05*rand(N + DiffEqFlux.paramlength(unn))
ps              = 0.01*rand(N + DiffEqFlux.paramlength(Hd))
ps[end-N+1:end] = 0.1*rand(N)
satu            = 1.0

include("dynamics.jl")
include("contactMap.jl")


sys = RimlessWheel(Float64, Δt, totalTimeStep)
cm  = ContactMap(Float64, sys, ps)
x0  = [0.2, cm.sys.γ, -2.0, 0.0]

# X, t = fulltimestep(cm, x0, Float64[])
# plots(X, t)

include("learnControl.jl")