using LinearAlgebra, Flux, DiffEqFlux
using PyPlot 

include("customPosterior.jl")
include("trainingHelpers.jl")
include("testHelpers.jl")

const tspan         = (0.0f0, 10.0f0) 
const batchsize     = 4
const Δt            = 0.01f0
const totalTimeStep = Int(floor(tspan[2]/Δt))
const binSize       = 4

nn = FastChain(FastDense(2, 4, elu),
                FastDense(4, 4))

inputLayer(x) = x

const nn_length = DiffEqFlux.paramlength(nn) 

function bin(x, θ) 
    return Flux.softmax(nn(inputLayer(x), θ))
end

function unstackParams(param)
    return param[1:nn_length], param[nn_length+1:nn_length+binSize]
end

function lossPerState(x)
    return 3.0f0*dot(x, x)
end

function computeLoss(x0, param::Vector{T}; totalTimeStep = totalTimeStep) where {T<:Real}
    
    ψ, θuk  = unstackParams(param)
    θuk_sig = berAct.(θuk)
    ltotal  = 0.0f0

    for xi in x0
        X  = Vector{Vector{T}}()
        pk = Vector(undef, totalTimeStep)

        push!(X, xi)

        for i in 1:totalTimeStep
            pk[i]   = bin(X[end], ψ)           #holds a softmax
            lk      = 0.0f0
            for k in 1:binSize

                lj = 0.0f0
                θk = [θuk_sig[k], 1.0f0-θuk_sig[k]]
                for j in 0:1 
                    x  = oneStep(X[end], j)
                    lj += θk[j+1]*lossPerState(x)       #maximize the changes of selecting lower loss
                end
                lk += pk[i][k]*lj

            end
            ltotal += lk/totalTimeStep

            ### select the next state
            c  = argmax(bin(X[end], ψ))
            uk = maxBernoulli(θuk[c])
            push!(X, oneStep(X[end], uk))
        end
    end

    return ltotal/length(x0)
end

function berAct(x)
    return Flux.sigmoid(x)
end

function maxBernoulli(θuk)
    if berAct(θuk) >= 0.5f0
        return 0.0f0 
    else 
        return 1.0f0 
    end
end

function trainEM()

    ψ           = randn(Float32, nn_length)
    θuk         = LogExpFunctions.logit.([0.6f0, 0.4f0, 0.4f0, 0.4f0])
    param       = vcat(ψ, θuk)
    opt         = Adam(0.001f0)
    counter     = 0
    minibatch   = 4

    for i in 1:5000
        uk  = maxBernoulli.(θuk)
        x0  = sampleInitialState(ψ, uk; totalTimeStep=2000, minibatch=minibatch)

        l1(θ)   = computeLoss(x0, θ; totalTimeStep = 1000)
        lg1     = ForwardDiff.gradient(l1, param)

        if counter > 10
            ψ, θuk  = unstackParams(param)
            uk      = maxBernoulli.(θuk)
            l       = testBayesian(x0[1], ψ, uk; totalTimeStep=1000)
            println("loss = ", l, " θuk = ", berAct.(θuk) ," grad = ", lg1[end])
            counter = 0
        end
        counter += 1
        Flux.update!(opt, param, lg1)
    end
end