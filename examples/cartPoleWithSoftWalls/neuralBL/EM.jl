using LinearAlgebra, Flux, DiffEqFlux, LogExpFunctions, ForwardDiff
using Distributions
using PyPlot 

using ContactLCP

include("../dynamics.jl")

sys  = CartPoleWithSoftWalls()
lcp  = ContactLCP.Lcp(Float32, sys)

include("trainingHelpers.jl")
include("testHelpers.jl")

const tspan         = (0.0f0, 4.0f0) 
const Δt            = 0.001f0
const totalTimeStep = Int(floor(tspan[2]/Δt))
const binSize       = 1

binNN = FastChain(FastDense(5, 6, elu),
                 FastDense(6, binSize))

const binNN_length = DiffEqFlux.paramlength(binNN) 

controlArray     = Array{Function}(undef, binSize);
controlNN_length = Vector{Int}(undef, binSize)
ps               = Vector{Vector}(undef, binSize)

for i in 1:binSize 
    controlArray[i] = FastChain(FastDense(5, 8, elu),
                            FastDense(8, 4, elu),
                            FastDense(4, 1))
    
    controlNN_length[i] = DiffEqFlux.paramlength(controlArray[i]) 

    ps[i]               = randn(Float32, controlNN_length[i])
end

function bin(x, θ::Vector{T}) where {T<:Real} 
    # return Flux.softmax(binNN(inputLayer(x), θ))
    return 1.0f0
end

function input(x::Vector{T}, θk::Vector{Vector}, i::Int) where {T<:Real}
    return controlArray[i](inputLayer(x), θk[i])
end

function unstackControlParam(param)

    θk = Vector{Vector}(undef, binSize)
    θk[1] = param[1:controlNN_length[1]] 

    for i in 2:binSize
        θk[i] = param[sum([controlNN_length[j] for j in 1:i-1])+1:sum([controlNN_length[j] for j in 1:i])] 
    end

    return θk
end

function unstackParams(param)
    return param[1:binNN_length], unstackControlParam(param[binNN_length+1:end])
end

function stackParams(ψ::Vector{T}, θk::Vector{Vector}) where {T<:Real}
    param = Vector{T}()
    param = vcat(param, ψ)
    for i in 1:binSize
        append!(param, θk[i])
    end
    return param
end

function lossPerState(x)
    x1, x2, x1dot, x2dot = x
    Q = diagm(0 => [0.0f0, 2.0f0, 2.0f0])

    return 1.0f0*dot([x1, x1dot, x2dot], 
                    Q*[x1, x1dot, x2dot]) + 4.0f0*(1.0f0-cos(x2))
end

function computeLoss(x0, param::Vector{T}; totalTimeStep = totalTimeStep) where {T<:Real}
    
    ψ, θk   = unstackParams(param)
    ltotal  = 0.0f0

    for xi in x0
        pk = Vector(undef, totalTimeStep)
        x2 = deepcopy(xi)

        for i in 1:totalTimeStep
            pk[i]  = bin(x2, ψ)           #holds a softmax
            lk     = 0.0f0
            for k in 1:binSize     
                u  = input(x2, θk, k)          
                x  = oneStep(x2, u)
                lk += pk[i][k]*lossPerState(x)       #maximize the changes of selecting lower loss
            end
            ltotal += lk/totalTimeStep

            ### select the next state
            c  = rand(Categorical(bin(x2, ψ)))
            # c  = argmax(bin(x2, ψ))
            x2 = oneStep(x2, input(x2, θk, c))
        end
    end

    return 3.0f0*ltotal/length(x0)
end

function trainEM()

    ψ           = randn(Float32, binNN_length)
    θk          = deepcopy(ps)
    param       = stackParams(ψ, θk)
    opt         = Adam(0.001f0)
    counter     = 0
    minibatch   = 2

    for i in 1:5000
        ψ, θk = unstackParams(param)
        x0    = sampleInitialState(ψ, θk; totalTimeStep=2000, minibatch=minibatch)

        l1(θ)   = computeLoss(x0, θ; totalTimeStep = 2000)
        lg1     = ForwardDiff.gradient(l1, param)

        if counter > 5
            ψ, θk  = unstackParams(param)
            testBayesian(x0[1], ψ, θk; totalTimeStep=4000)
            println("loss = ", l1(param))
            counter = 0
        end
        counter += 1
        Flux.update!(opt, param, lg1)
    end
end

function checkGradient(param)
    x0  = [[0.0f0, pi, 1.0f0, 0.5f0]]

    l1(θ)  = computeLoss(x0, θ; totalTimeStep = 1000)
    lg1    = ForwardDiff.gradient(l1, param)

    δ = 0.001
    fd = zeros(length(param))

    for ind in 1:length(param)
        param2 = deepcopy(param)
        param2[ind] += δ

        fd[ind] = (l1(param2) - l1(param))/δ
    end

    percError = abs.(lg1 - fd) ./ fd

    println("error = ", percError)
end