using DiffEqFlux, MLBasedESC, Flux
using LinearAlgebra, Random
using CUDA

rng = Random.default_rng()
Random.seed!(rng, 0)

include("dynamics.jl")
include("learnControl.jl")

const m1  = 2.32f0     
const m2  = 4.194f0 
const I1  = 0.0784160f0
const I2  = 0.0380256f0
const mt  = m1 + m2
const l1  = 0.26f0
const l2  = 0.05f0
const g   = 9.81f0
const α   = 360.0f0/10.0f0/2.0f0 * pi/180.0f0
const γ   = 0.0f0*pi/180.0f0
const ϵn  = 0.0f0
const Δt = 0.001f0
const totalTimeStep = 5.0f0/Δt
const total_contact_num = 1
const gThreshold = 0.001f0


Hd              = FastChain(FastDense(6, 8, elu), FastDense(8, 5, elu), FastDense(5, 1))
N               = 6
npbc            = MLBasedESC.NeuralPBC(N, Hd)

ps              = CuArray(0.5f0*CUDA.randn(Float32, N + DiffEqFlux.paramlength(Hd)))
ps[end-N+1:end] = 0.1f0*CUDA.rand(Float32, N)
satu            = 2.0f0
x0              = [0.2f0, 0.0f0, -2.0f0, 0.0f0]