using Turing, Distributions, OrdinaryDiffEq, DistributionsAD
using LinearAlgebra, Flux, DiffResults, DiffEqFlux
using ProgressMeter, AdvancedVI, Distributions, LogExpFunctions
using Turing: Variational, ForwardDiff
using PyPlot 

include("customPosterior.jl")
include("trainingHelpers.jl")
include("testHelpers.jl")

const tspan         = (0.0f0, 10.0f0) 
const batchsize     = 4
const Δt            = 0.01f0
const totalTimeStep = Int(floor(tspan[2]/Δt))
const binSize       = 4

nn = FastChain(FastDense(2, 8, elu),
                FastDense(8, 4))

# inputLayer(x) = [cos(x[1]), sin(x[1]), cos(x[2]), sin(x[2])]
inputLayer(x) = x

const nn_length = DiffEqFlux.paramlength(nn) 

function unstackDistParam(par)
    # return par[1:nn_length], par[nn_length+1:2nn_length], par[2nn_length+1:2nn_length+binSize]
    return par[1:binSize]
end

function unstackSamples(par)
    # return par[1:nn_length], par[nn_length+1:nn_length+binSize]
    return par[1:binSize]
end

function trajectory(x0, par; totalTimeStep = totalTimeStep)
    θ      = rand(q(par))
    ψ, uk  = unstackSamples(θ)

    trajectory(x0, ψ, uk; totalTimeStep = totalTimeStep)
    return X
end

function sampleInitialState(par::Vector{T}; totalTimeStep = totalTimeStep, minibatch=4) where {T<:Real}
    
    x0 = [rand(-4.0:0.001:4.0), rand(-4.0:0.001:4.0)]
    X  = trajectory(x0, par)
    
    samples = Vector{Vector{T}}()
    for i in 1:minibatch
        if rand() > 0.1
            push!(samples, rand(X))
        else
            push!(samples, [rand(-0.1:0.001:0.1), rand(-0.1:0.001:0.1)])
        end

    end

    return samples
    # return x0
end

function bin(x, θ) 
    return Flux.softmax(nn(inputLayer(x), θ))
end

function loss(traj)    
    # set_distance_loss(traj)
    accumulatedLoss(traj)
end

function binLoss(x0, ψ::Vector{T}, uk; totalTimeStep = totalTimeStep) where {T<:Real}
    l = 0.0f0
    for xi in x0
        l  += computeLoss(xi, ψ, uk; totalTimeStep = totalTimeStep)
    end
    return l/length(x0)
end

function q(par::Vector{T}) where {T<:Real}
    θk = unstackDistParam(par)
    return PosteriorDistribution(θk)
    # return MvNormal(θk, 0.001)
end

@model function stateBinModel(x0, par::Vector{T}; binSize = 4, totalTimeStep = totalTimeStep) where {T<:Real}

    #define ψ that parameterizes the function that takes states and returns categorical probabilities
    mψ, σψ, θk = unstackDistParam(par)

    ψ   = Vector{T}(undef, nn_length)
    ψ   ~ MvNormal(mψ, Flux.softplus.(σψ))   #currently deterministic
    
    uk  = Vector{T}(undef, binSize)
    for i in 1:binSize
        uk[i] ~ Bernoulli(berAct(θk[i]))     #for one whole trajectory, each bin only gets one sample of control. Do not resample controller for each bin for every state 
    end

    for xi in x0
        X,_   = integrate(xi, ψ, uk; totalTimeStep = totalTimeStep)
        l     = loss(X)
        0.0f0 ~ Normal(l, 2.0f0)
    end

end

@model function stateBinModel(x0, par::Vector{T}, ψ; binSize = 4, totalTimeStep = totalTimeStep) where {T<:Real}

    #define ψ that parameterizes the function that takes states and returns categorical probabilities
    # θk = unstackDistParam(par)
    uk  = Vector{T}(undef, binSize)
    # uk ~ arraydist(Bernoulli.(berAct.(θk)))     #for one whole trajectory, each bin only gets one sample of control. Do not resample controller for each bin for every state 

    w   ~ Dirichlet(binSize, 1.0)
    uk  ~ arraydist(Bernoulli.(w))

    l = 0.0f0
    for xi in x0
        l     +=  computeLoss(xi, ψ, uk; totalTimeStep = totalTimeStep)
        0.0f0 ~ Normal(l, 0.5f0)
    end
end

function trainDeterPsi()
    
    ψ           = randn(Float32, nn_length)
    θk          = LogExpFunctions.logit.([0.8f0, 0.8f0, 0.3f0, 0.3f0])   
    vo          = Variational.ELBO()
    elbo_num    = 5
    alg         = ADVI(elbo_num, 100)

    optimizer   = Variational.DecayedADAGrad(0.001f0)
    opt         = Adam(0.001f0)
    par         = deepcopy(θk)
    diff_results = DiffResults.GradientResult(par)

    counter     = 0
    minibatch   = 4
    fullChain = []

    actions = []
    θuk = []
    for i in 1:binSize 
        push!(actions, FastChain(FastDense(2, 4, elu),
                        FastDense(4, 2)))
        push!(θuk, DiffEqFlux.paramlength(actions[i]))
    end

    for i in 1:5000

        x0   = sampleInitialState(ψ, uk; totalTimeStep=2000, minibatch=minibatch)

        l1(param) = binLoss(x0, param, uk; totalTimeStep = 1000)
        lg1 = ForwardDiff.gradient(l1, ψ)
        Flux.update!(opt, ψ, lg1)

        l2(param) = controlLoss(x0, param, θuk; totalTimeStep = 1000)
        lg2 = ForwardDiff.gradient(l2, θuk)
        Flux.update!(opt, θuk, lg2)
        # model   = stateBinModel(x0, par, ψ; totalTimeStep=1000)
        # AdvancedVI.grad!(vo, alg, q, model, par, diff_results, elbo_num)
        # ∇       = DiffResults.gradient(diff_results)

        if counter > 10
            l   = testBayesian(x0[1], ψ, uk ; totalTimeStep=1000)
            println("loss = ", l, " uk = ",uk)
            counter = 0
        end
        # Δ       = Flux.Optimise.apply!(optimizer, par, ∇)
        # @. par  = par - Δ
      
        # sampler = Gibbs(PG(10, :uk, :w))
        # nsamples = 100
        # chains = sample(model, IS(), nsamples; progress=false);
        # uk  = [mode(chains[["uk[$i]"]].value.data) for i in 1:binSize]
        # w   = [mean(chains[["w[$i]"]].value.data) for i in 1:binSize]
        # push!([uk, w])
        counter += 1
    end
end

function trainBayesianStateBinModel()
    
    mψ          = randn(Float32, nn_length)
    σψ          = LogExpFunctions.invsoftplus.(0.1 .* rand(Float32, nn_length))
    θk          = LogExpFunctions.logit.([0.8f0, 0.8f0, 0.3f0, 0.3f0])   
    vo          = Variational.ELBO()
    elbo_num    = 5
    alg         = ADVI(elbo_num, 100)

    optimizer   = Variational.DecayedADAGrad(0.001f0)
    par         = vcat(mψ, σψ, θk)
    diff_results = DiffResults.GradientResult(par)

    counter     = 0
    minibatch   = 4

    for i in 1:5000

        x0      = sampleInitialState(par; totalTimeStep=2000, minibatch=minibatch)

        model   = stateBinModel(x0, par; totalTimeStep=1000)
        AdvancedVI.grad!(vo, alg, q, model, par, diff_results, elbo_num)
        ∇       = DiffResults.gradient(diff_results)

        if counter > 10
            l   = testBayesian(x0[1], par; totalTimeStep=1000)
            println("loss = ", l, " θk = ", berAct.(unstackDistParam(par)[3]), " grad = ",  ∇[end])
            counter = 0
        end
        Δ       = Flux.Optimise.apply!(optimizer, par, ∇)
        @. par  = par - Δ
      
        counter += 1
    end
end

function trainRapidSwitching()  

    μθ       = rand(Float32, nn_length)
    σθ       = 0.01f0*ones(Float32, nn_length)

    vo          = Variational.ELBO()
    elbo_num    = 5
    alg         = ADVI(elbo_num, 100)

    optimizer   = Variational.DecayedADAGrad(0.001f0)
    θdistParam   = vcat(μθ, σθ)
    diff_results = DiffResults.GradientResult(param_cur)

    counter     = 0
    minibatch   = 1

    for i in 1:5000

        m_cur, standev_cur = unstack(θdistParam)
        x0 = [0.5f0, 0.5f0]
        # x0 = sampleInitialState(θdistParam; totalTimeStep=1000, minibatch=minibatch)

        model   = rapidSwitchingModel(x0, m_cur, standev_cur; totalTimeStep=500)
        AdvancedVI.grad!(vo, alg, q, model, θdistParam, diff_results, elbo_num)
        ∇       = DiffResults.gradient(diff_results)

        if counter > 10
            l   = testTrajectory(x0, θdistParam; totalTimeStep=500)
            println("loss = ", l, " grad = ",  ∇[1])
            counter = 0
        end
        Δ       = Flux.Optimise.apply!(optimizer, θdistParam, ∇)
        @. θdistParam = θdistParam - Δ
      
        counter += 1
    end
end

