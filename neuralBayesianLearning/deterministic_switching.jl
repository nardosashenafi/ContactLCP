using LinearAlgebra, Flux, DiffEqFlux
using PyPlot 

include("customPosterior.jl")

const tspan         = (0.0f0, 10.0f0) 
const batchsize     = 4
const Δt            = 0.01f0
const totalTimeStep = Int(floor(tspan[2]/Δt))
const binSize       = 4

nn = FastChain(FastDense(2, 8, elu),
                FastDense(8, 4))


function unstackDistParam(par)
    return par[1:nn_length], par[nn_length+1:2nn_length], par[2nn_length+1:2nn_length+binSize]
end

function unstackSamples(par)
    return par[1:nn_length], par[nn_length+1:nn_length+binSize]
end

function trajectory(x0, par; totalTimeStep = totalTimeStep)
    θ      = rand(q(par))
    ψ, uk  = unstackSamples(θ)
    
    X,_    = integrate(x0, ψ, uk; totalTimeStep = totalTimeStep)

    return X
end

function bin(x, θ) 
    return Flux.softmax(nn(inputLayer(x), θ))
end


function trainDeterministicModel()
    mψ          = randn(Float32, nn_length)

    for i in 1:5000

    end
    
end
