using Turing, Distributions, OrdinaryDiffEq, DistributionsAD
using Turing: Variational
using ContactLCP

include("../dynamics.jl")
include("trainingHelpers.jl")
include("testHelpers.jl")

const tspan         = (0.0f0, 10.0f0) 
const batchsize     = 4
const Δt            = 0.01f0
const totalTimeStep = Int(floor(tspan[2]/Δt))
const binSize       = 4

binNN = FastChain(FastDense(5, 6, elu),
                 FastDense(6, 3))

const binNN_length = DiffEqFlux.paramlength(binNN) 

controlArray     = Array{Function}(undef, binSize);
controlNN_length = Vector{Int}(undef, binSize)
ps               = Vector{Vector}(undef, binSize)

for i in 1:binSize 
    controlArray[i] = FastChain(FastDense(5, 6, elu),
                        FastDense(6, 1))
    
    controlNN_length[i] = DiffEqFlux.paramlength(controlArray[i]) 

    ps[i]               = randn(Float32, controllNN_length[i])
end

function bin(x, θ) 
    return Flux.softmax(binNN(inputLayer(x), θ))
end

function unstackDistParam(par)
    return par[1:binSize]
end

function unstackSamples(par)
    return par[1:binSize]
end

function lossPerState(x)
    return 3.0f0*dot(x, x)
end

function q(par::Vector{T}) where {T<:Real}
    θk = unstackDistParam(par)
    return PosteriorDistribution(θk)
    # return arraydist(Normal.(berAct.(θk)))
end

function computeLoss(x0, ψ, uk::Vector{T}; totalTimeStep = totalTimeStep) where {T<:Real}
    
    ltotal  = 0.0f0

    for xi in x0
        X  = Vector{Vector{T}}()
        pk = Vector(undef, totalTimeStep)

        push!(X, xi)

        for i in 1:totalTimeStep
            pk[i]   = bin(X[end], ψ)           #holds a softmax
            lk      = 0.0f0
            for k in 1:binSize
                x  = oneStep(X[end], uk[k])
                lk += pk[i][k]*lossPerState(x)
            end
            ltotal += lk/totalTimeStep

            ### select the next state
            c  = rand(Categorical(bin(X[end], ψ)))
            push!(X, oneStep(X[end], uk[c]))
        end
    end

    return ltotal/length(x0)
end

@model function stateBinModel(x0, par::Vector{T}, ψ; binSize = 4, totalTimeStep = totalTimeStep) where {T<:Real}

    θk    = unstackDistParam(par)
    uk    = Vector{T}(undef, binSize)
    for i in 1:binSize
        uk[i] ~ Bernoulli(berAct(θk[i]))     #for one whole trajectory, each bin only gets one sample of control. Do not resample controller for each bin for every state 
    end
    l     =  computeLoss(x0, ψ, uk; totalTimeStep = totalTimeStep)
    0.0f0 ~ Normal(l, 0.5f0)

end

##############################
#Currently, VI is not working because Advanced.grad! is not taking gradient of the 
#bernoulli sample wrt its distribution parameters θk accurately. Needs proper definition
#of a bijector that transforms the bernoulli to a continuous function such as 
#Gumbel-softmax

function trainVI()

    ψ           = randn(Float32, nn_length)
    θk          = LogExpFunctions.logit.([0.8f0, 0.3f0, 0.3f0, 0.3f0])   
    vo          = Variational.ELBO()
    elbo_num    = 5
    alg         = ADVI(elbo_num, 100)

    optimizer   = Variational.DecayedADAGrad(0.001f0)
    opt         = Adam(0.001f0)
    par         = deepcopy(θk)
    diff_results = DiffResults.GradientResult(par)

    counter     = 0
    minibatch   = 4

    for i in 1:5000
        uk   = rand(q(par))
        x0   = sampleInitialState(ψ, uk; totalTimeStep=2000, minibatch=minibatch)

        l1(param)   = computeLoss(x0, param, uk; totalTimeStep = 1000)
        lg1         = ForwardDiff.gradient(l1, ψ)
        Flux.update!(opt, ψ, lg1)

        model   = stateBinModel(x0, par, ψ; totalTimeStep=1000)
        AdvancedVI.grad!(vo, alg, q, model, par, diff_results, elbo_num)
        ∇       = DiffResults.gradient(diff_results)

        if counter > 10
            testBayesian(x0[1], ψ, uk ; totalTimeStep=1000)
            println("loss = ", l1(ψ), " par = ", berAct.(par))
            counter = 0
        end
        Δ       = Flux.Optimise.apply!(optimizer, par, ∇)
        @. par  = par - Δ

        counter += 1
    end

end
