
using DiffEqFlux
using MLBasedESC
using Statistics
using ProgressMeter
using MeshCat, GeometryBasics, CoordinateTransformations, ColorTypes, Blink, Rotations

include("bayesianTrainingHelpers.jl")

Δt = 0.001f0; totalTimeStep = 1500

x0             = Float32.(initialState(pi, -1.0f0, 0.0f0, 0.0f0))
param_expert   = Float32[30.0, 5.0]

Hd              = FastChain(FastDense(6, 8, elu), 
                  FastDense(8, 5, elu),
                  FastDense(5, 1))
const N         = 6
npbc            = MLBasedESC.NeuralPBC(N, Hd)
const satu     = 2.0f0
