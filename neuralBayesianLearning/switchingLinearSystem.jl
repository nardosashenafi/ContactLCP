using Turing, Distributions, OrdinaryDiffEq, DistributionsAD
using LinearAlgebra, Flux, DiffResults, DiffEqFlux
using ProgressMeter, AdvancedVI, Distributions
using Turing: Variational
using PyPlot 

include("customPosterior.jl")

const tspan         = (0.0f0, 10.0f0) 
const batchsize     = 4
const Δt            = 0.01f0
const totalTimeStep = Int(floor(tspan[2]/Δt))
const binSize       = 4

nn = FastChain(FastDense(2, 5, elu),
                FastDense(5, 5, elu),
                FastDense(5, 4))

const nn_length = DiffEqFlux.paramlength(nn) 

function oneStep(x, ρ::T) where {T<:Real}
    A1 = [0.0f0 -1.0f0; 2.0f0 0.0f0]
    A2 = [0.0f0 -2.0f0; 1.0f0 0.0f0]

    dx = (1.0f0 - ρ)*A1*x + ρ*A2*x
    return x + dx*Δt
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

function trajectory(x0, par)
    θ      = rand(q(par))
    ψ, uk  = unstackSamples(θ)
    
    pk = Vector(undef, totalTimeStep)
    c = Vector{Int}(undef, totalTimeStep)

    X = Vector{Vector{T}}()
    push!(X, x0)

    for i in 1:totalTimeStep
        pk[i] = bin(X[end], ψ)           #holds a softmax
        c[i]  = argmax(pk[i])
        x     = oneStep(X[end], uk[c[i]])
        push!(X, x)
    end
    return X
end

function loss(traj)

    extract = traj[end-100:end]
    l = 5.0f0/length(extract)*sum(map(x -> dot(x, x), extract))
    # @assert l > 1.0f0 "loss is not normalized" 
end

function unstackDistParam(par)
    return par[1:nn_length], par[nn_length+1:nn_length+binSize]
end

function unstackSamples(par)
    return par[1:nn_length], par[nn_length+1:nn_length+binSize]
end

function q(par::Vector{T}) where {T<:Real}
    # m, standev = unstack(p)
    mψ, θk = unstackDistParam(par)
    return PosteriorDistribution(mψ, 0.00001f0*ones(T, nn_length), θk)
end

function bin(x, θ) 
    # return θ[1]*x[1]^2 + θ[2]*x[2]^2 + θ[3]*x[1]*x[2]
    return softmax(nn(x, θ))
end

@model function stateBinModel(x0::Vector{Vector{T}}, par; binSize = 4, totalTimeStep = totalTimeStep) where {T<:Real}

    #define ψ that parameterizes the function that takes states and returns categorical probabilities
    mψ, θk = unstackDistParam(par)

    ψ = Vector{T}(undef, nn_length)
    ψ ~ MvNormal(mψ, 1e-5)   #currently deterministic
    
    uk = Vector{T}(undef, binSize)
    for i in 1:binSize
        uk[i] ~ Bernoulli(Flux.sigmoid(θk[i]))     #for one whole trajectory, each bin only gets one sample of control. Do not resample controller for each bin for every state 
    end
    pk = Vector(undef, totalTimeStep)
    c = Vector{Int}(undef, totalTimeStep)

    for xi in x0
        X = Vector{Vector{T}}()
        push!(X, xi)

        for i in 1:totalTimeStep
            pk[i] = bin(X[end], ψ)           #holds a softmax
            c[i]  = argmax(pk[i])
            x     = oneStep(X[end], uk[c[i]])
            push!(X, x)
        end

        l     = loss(X)
        0.0f0 ~ Normal(l, 0.1f0)
    end

end

function sampleInitialState(par::Vector{T}; totalTimeStep = totalTimeStep, minibatch=4) where {T<:Real}
    
    x0 = [rand(-2.0:0.01:2.0), rand(-2.0:0.01:2.0)]
    X  = trajectory(x0, par)
    
    samples = Vector{Vector{T}}()
    for i in 1:minibatch
        if rand() > 0.2
            push!(samples, rand(X))
        else
            push!(samples, [rand(-0.1:0.001:0.1), rand(-0.1:0.001:0.1)])
        end

    end

    return samples
    # return x0
end

function trainBayesianStateBinModel()
    
    mψ          = rand(Float32, nn_length)
    θk          = [0.8f0, 0.8f0, 0.2f0, 0.2f0]      
    vo          = Variational.ELBO()
    elbo_num    = 7
    alg         = ADVI(elbo_num, 100)

    optimizer   = Variational.DecayedADAGrad(0.005f0)
    par         = vcat(mψ, θk)
    diff_results = DiffResults.GradientResult(par)

    counter     = 0
    minibatch   = 2

    for i in 1:5000

        x0 = sampleInitialState(par; totalTimeStep=500, minibatch=minibatch)

        model   = stateBinModel(x0, par; totalTimeStep=500)
        AdvancedVI.grad!(vo, alg, q, model, par, diff_results, elbo_num)
        ∇       = DiffResults.gradient(diff_results)

        if counter > 20
            l, _   = testBayesian(x0[1], par; totalTimeStep=1000)
            println("loss = ", l, " θk = ", unstackDistParam(par)[2], " grad = ",  ∇[1])
            counter = 0
        end
        Δ       = Flux.Optimise.apply!(optimizer, par, ∇)
        @. par  = par - Δ
      
        counter += 1
    end
end


function testBayesian(xi, par; totalTimeStep = totalTimeStep)
    θ     = rand(q(par))
    ψ, uk = unstackSamples(θ)
    
    pk = Vector(undef, totalTimeStep)
    c = Vector{Int}(undef, totalTimeStep)

    X = Vector{Vector{T}}()
    U = Vector{T}()
    push!(X, xi)

    for i in 1:totalTimeStep
        pk[i] = bin(X[end], ψ)           #holds a softmax
        c[i]  = argmax(pk[i])
        x     = oneStep(X[end], uk[c[i]])
        push!(X, deepcopy(x))
        push!(U, deepcopy(uk[c[i]]))
    end
    clf()
    plt.subplot(2, 2, 1)
    plot(getindex.(X, 1), getindex.(X, 2))
    scatter(X[end][1], X[end][2], marker="o")
    l = loss(X)
    ylabel("x2")
    xlabel("x1")
    legend(["loss " * string(l)])

    plt.subplot(2, 2, 2)
    plot(getindex.(X[1:end-1], 1), U)
    ylabel("uk")

    plotPartition(ψ, uk, X)
    return l, U
end

function plotPartition(ψ, uk, X)
    width = 20
    X1 = range(-2.0, 2.0, length=width)
    X2 = range(-2.0, 2.0, length=width)
    
    u = Matrix{Int}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            pk      = bin(vcat(X1[j], X2[width-i+1]), ψ)
            c[i, j] = argmax(pk)
            u[i, j] = uk[c[i, j]]
        end
    end

    plt.subplot(2, 2, 3)
    plot(getindex.(X, 1), getindex.(X, 2))
    scatter(X[end][1], X[end][2], marker="o")
    l = loss(X)
    ylabel("x2")
    xlabel("x1")
    legend(["loss " * string(l)])

    imshow(u, extent = [X1[1], X1[end], X2[1], X2[end]])
    ylabel("Control")

    plt.subplot(2, 2, 4)
    imshow(c, extent = [X1[1], X1[end], X2[1], X2[end]])
    ylabel("Partitions")
end


function checkGradient()
    mψ          = rand(Float32, nn_length)
    θk          = [0.5f0, 0.5f0, 0.5f0, 0.5f0]      
    vo          = Variational.ELBO()
    elbo_num    = 10
    alg         = ADVI(elbo_num, 100)

    par         = vcat(mψ, θk)
    diff_results = DiffResults.GradientResult(par)

    minibatch   = 2

    x0 = sampleInitialState(par; totalTimeStep=500, minibatch=minibatch)

    model   = stateBinModel(x0, par; totalTimeStep=500)
    AdvancedVI.grad!(vo, alg, q, model, par, diff_results, elbo_num)
    ∇       = DiffResults.gradient(diff_results)

    par2 = deepcopy(par)
    par2[1] += 0.001
    
    

end

function evaluateTraining(chains, x0)
    b = [chains[["b[$i]"]].value.data[end] for i in 1:totalTimeStep]
    u = [chains[["u[$i]"]].value.data[end] for i in 1:binSize]
    t = [chains[["transition[$i]"]].value.data[end] for i in 1:binSize]
    θ = [chains[["θ[$i]"]].value.data[end] for i in 1:nn_length]

    X = trajectory(x0, p)

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

