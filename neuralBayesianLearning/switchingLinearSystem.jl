using Turing, Distributions, OrdinaryDiffEq, DistributionsAD
using LinearAlgebra, Flux, DiffResults, DiffEqFlux
using ProgressMeter, AdvancedVI
using Turing: Variational

const tspan         = (0.0f0, 10.0f0) 
const batchsize     = 4
const Δt            = 0.05f0
const totalTimeStep = Int(floor(tspan[2]/Δt))

nn = FastChain(FastDense(2, 5, elu),
        FastDense(5, 2))

const nn_length = DiffEqFlux.paramlength(nn) 

function oneStep(x, ρ::T) where {T<:Real}
    A1 = [0.0f0 -1.0f0; 2.0f0 0.0f0]
    A2 = [0.0f0 -2.0f0; 1.0f0 0.0f0]

    dx = ρ*A1*x + (1.0f0-ρ)*A2*x
    return x + dx*Δt
end

function trajectory(x0, θ)
    X = Vector{Vector{T}}()
    push!(X, x0)

    for i in 1:totalTimeStep-1
        p   = categorize(X[end], θ)
        x    = oneStep(X[end], p)
        push!(X, x)
    end

    return X
end


function loss(traj)
    l = 1.0f0/length(traj)*sum(map(x -> dot(x, x), traj))
    # @assert l > 1.0f0 "loss is not normalized" 
end

function f(x, θ) 
    # return θ[1]*x[1]^2 + θ[2]*x[2]^2 + θ[3]*x[1]*x[2]
    return nn(x, θ)
end

function unstack(p)
    return p[1:nn_length], p[nn_length+1:nn_length*2]
end

function q(p)
    m, standev = unstack(p)
    return MvNormal(m, Flux.softplus.(standev))
end

function desiredTrajectory(x0)
    X     = []
    push!(X, x0)

    for i in 1:totalTimeStep-1
        if prod(X[end]) <= 0
            p = 1
        else
            p = 0
        end
        x    = oneStep(X[end], p)
        push!(X, x)
    end
    clf()
    plot(getindex.(X, 1), getindex.(X, 2))
    l = loss(X)
    return l
end

function testTrajectory(x0, par; totalTimeStep = totalTimeStep)
    theta = rand(q(par))
    X     = []
    push!(X, x0)

    for i in 1:totalTimeStep-1
        p    = categorize(X[end], theta)
        x    = oneStep(X[end], p)
        push!(X, x)
    end
    clf()
    plot(getindex.(X, 1), getindex.(X, 2))
    l = loss(X)
    return l
end

function categorize(x, θ)
    # return rand(Bernoulli(f(X[end], θ)))
    argmax(softmax(f(x, θ))) - 1
    # softmax(f(x, θ))
end

@model function switchingModel(x0::Vector{Vector{T}}, μθ, σθ; totalTimeStep = totalTimeStep) where {T<:Real}

    #define θ that parameterizes the function that takes states and returns categorical probabilities
    θ = Vector{T}(undef, nn_length)
    θ ~ MvNormal(μθ, Flux.softplus.(σθ))

    p = Vector{T}(undef, totalTimeStep-1)

    for xi in x0
        X = Vector{Vector{T}}()
        push!(X, xi)

        for i in 1:totalTimeStep-1
            p[i] = categorize(X[end], θ)
            # p[i] ~ Categorical(softmax(f(X[end], θ)))
            x    = oneStep(X[end], p[i]-1)
            push!(X, x)
        end

        l = loss(X)
        # println("loss = ", l)
        0.0f0 ~ Normal(l, 0.1f0)
    end

end

function sampleInitialState(θ::Vector{T}; totalTimeStep = totalTimeStep, minibatch=4)
    
    x0 = [rand(-1.5:0.1:1.5), rand(1.5:0.1:1.5)]
    X  = trajectory(x0, θ)
    
    samples = Vector{Vector{T}}()
    for i in 1:minibatch
        push!(samples, rand(X))
    end

    return samples
end

function train()  
    # x0 = [0.5f0, 0.5f0] 
    # sampler = HMC(0.05, 100)
    # sampler = MH()
    # nsamples = 1000
    # nchains = 3
    # chains = sample(model, sampler, nsamples)

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
        # x0 = [rand(-1.5:0.1:1.5), rand(1.5:0.1:1.5)]
        x0 = sampleInitialState(θdistParam; totalTimeStep=1000, minibatch=minibatch)

        model   = switchingModel(x0, m_cur, standev_cur; totalTimeStep=Int(floor(3.0/Δt)))
        AdvancedVI.grad!(vo, alg, q, model, θdistParam, diff_results, elbo_num)
        ∇       = DiffResults.gradient(diff_results)

        if counter > 10
            l   = testTrajectory(x0[1], θdistParam; totalTimeStep=1000)
            println("loss = ", l, " grad = ",  ∇[1])
            counter = 0
        end
        Δ       = Flux.Optimise.apply!(optimizer, θdistParam, ∇)
        @. θdistParam = θdistParam - Δ
      
        counter += 1
    end
end

function evaluateTraining(chains, x0)
    theta = [chains[["θ[$i]"]].value.data[end] for i in 1:3]
    p = [chains[["p[$i]"]].value.data[end] for i in 1:totalTimeStep-1]

    X = trajectory(x0, p)

end
