
using DiffEqFlux, LogExpFunctions
using MLBasedESC
using Statistics
using ProgressMeter, DSP, BSON
using MeshCat, GeometryBasics, CoordinateTransformations, ColorTypes, Blink, Rotations
using Turing, AdvancedVI, Distributions 
using Turing: Variational 


include("bayesianTrainingHelpers.jl")

x0             = Float32.(initialState(pi, -1.0f0, 0.0f0, 0.0f0))
param_expert   = Float32[30.0, 5.0]

Hd              = FastChain(FastDense(6, 12, elu), 
                  FastDense(12, 7, elu),
                  FastDense(7, 1))
const N         = 6
const paramNum  = DiffEqFlux.paramlength(Hd)+N
npbc            = MLBasedESC.NeuralPBC(N, Hd)
const satu      = 2.0f0
responsetype    = DSP.Lowpass(0.5)
designmethod    = DSP.Butterworth(4)


@model function fitLoss(X0, R, μ_prior::AbstractArray{T}, σ_prior::AbstractArray{T}; totalTimeStep = 2000) where {T<:Real}

    w ~ arraydist(map((m, s) -> Distributions.Normal(m, LogExpFunctions.softplus(s)), μ_prior, σ_prior))

    l11 = 0.0f0
    for i in eachindex(X0)
        X , obstacles  = trajectory(X0[i], R[i], w; totalTimeStep=totalTimeStep)  
        l11 += hipSpeedLoss(X, obstacles)
    end
    
    0.0 ~ Distributions.Normal(l11, 0.001f0)
end

function getq(param)
    μ_param, σ_param = unstackParams(param)
    return MvNormal(μ_param, LogExpFunctions.softplus.(σ_param))
end

function controlToHipSpeed()

    fig1 = figure()
    fig2, ax2 = subplots()
    Turing.setrdcache(true)
    Turing.setadbackend(:forwarddiff)

    num_steps   = 5000
    elbo_num    = 3
    minibatchsize  = 2

    param = Float32.(vcat(0.3f0*randn(DiffEqFlux.paramlength(Hd)), 
                    0.1f0*rand(N), 
                    invsoftplus.(0.1f0*rand(DiffEqFlux.paramlength(Hd))),
                    invsoftplus.(0.05f0*rand(N))))

    optimizer   = Variational.DecayedADAGrad(0.01)

    vo          = Variational.ELBO()
    alg         = ADVI(elbo_num, num_steps)
    converged   = false

    prog        = ProgressMeter.Progress(num_steps, 1)
    diff_results = DiffResults.GradientResult(param) 
    elbo_data   = Vector{Float32}()
    X0, R       = sampleInitialStates(param, minibatchsize; totalTime=7000)
    μ_param, σ_param = unstackParams(param)
    model       = fitLoss(X0, R, μ_param, σ_param; totalTimeStep=4000)

    @showprogress for i in length(elbo_data)+1:5000

        X0, R       = sampleInitialStates(param, minibatchsize; totalTime=7000)
        println("X0 = ", X0)
        println("R = ", R)
        AdvancedVI.grad!(vo, alg, getq, model, param, diff_results, elbo_num)

        ∇ = DiffResults.gradient(diff_results)
        Δ = Flux.Optimise.apply!(optimizer, param, ∇)
        if any(isnan.(Δ))
            println("NaN occurred")
            continue
        end
        @. param = param - Δ
        ###################################################

        push!(elbo_data, vo(alg, getq(param), model, elbo_num))

        if mod(i, 5) == 0.0 
            converged = hasconverged(X0[1], R[1], param, elbo_data,i) 
            μ_param, σ_param = unstackParams(param)
            model       = fitLoss(X0, R, μ_param, σ_param; totalTimeStep=4000)
        end

        ProgressMeter.next!(prog)
    end
end
