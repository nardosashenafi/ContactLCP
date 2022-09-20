using LaTeXStrings, PyPlot
# using JuMP, PATHSolver, Mosek, MosekTools, DiffOpt
using LinearAlgebra
# import ChainRulesCore
using ForwardDiff, DiffEqFlux, ReverseDiff
using PyCall
@pyimport matplotlib.animation as anim
using Flux, Printf
using Interpolations
using Statistics
using MLBasedESC

Δt              = 0.001 
totalTimeStep   = 500

unn             = FastChain(FastDense(6, 8, elu), FastDense(8, 1))
Hd              = FastChain(FastDense(6, 8, elu), FastDense(8, 5, elu), FastDense(5, 1))
N               = 6
npbc            = MLBasedESC.NeuralPBC(N, Hd)

# println("NN parameterized control")
# ps              = 0.005*randn(N + DiffEqFlux.paramlength(unn))
# satu            = 1.5

println("NeuralPBC")
ps              = 0.5*randn(N + DiffEqFlux.paramlength(Hd))
ps[end-N+1:end] = 0.1*rand(N)
satu            = 2.0

include("dynamics.jl")
include("contactMap.jl")


sys = RimlessWheel(Float64, Δt, totalTimeStep)
cm  = ContactMap(Float64, sys, ps)
x0  = [0.2, cm.sys.γ, -2.0, 0.0]

# X, t = fulltimestep(cm, x0, Float64[])
# plots(X, t)

include("learnControl.jl")