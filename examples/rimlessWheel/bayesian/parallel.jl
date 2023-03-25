
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

Hd              = FastChain(FastDense(6, 8, elu), 
                  FastDense(8, 5, elu),
                  FastDense(5, 1))
const N         = 6
const paramNum  = DiffEqFlux.paramlength(Hd)+N
npbc            = MLBasedESC.NeuralPBC(N, Hd)
const satu      = 1.5f0
responsetype    = DSP.Lowpass(0.5)
designmethod    = DSP.Butterworth(4)

@model function fitLoss(x0, r, param::AbstractArray{T}; totalTimeStep = 2000) where {T<:Real}

    μ_prior, σ_prior = unstackParams(param)
    w ~ arraydist(map((m, s) -> Distributions.Normal(m, LogExpFunctions.softplus(s)), μ_prior, σ_prior))

    X, obstacles  = trajectory(x0, r, w; totalTimeStep=totalTimeStep)   
    l11     = hipSpeedLoss(X, obstacles)
    
    0.0 ~ Distributions.Normal(l11, 0.0001f0)
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

    optimizer   = Variational.DecayedADAGrad(0.005)

    vo          = Variational.ELBO()
    alg         = ADVI(elbo_num, num_steps)
    converged   = false

    prog        = ProgressMeter.Progress(num_steps, 1)
    diff_results = [DiffResults.GradientResult(param) for i in 1:minibatchsize]
    elbo_data   = Vector{Float32}()

    @showprogress for i in 1:5000

        X0, R = sampleInitialStates(param, minibatchsize; totalTime=7000)
        println("X0 = ", X0)
        println("R = ", R)
        Threads.@threads for i in eachindex(X0)
            model   = fitLoss(X0[i], R[i], param; totalTimeStep=4000)
            AdvancedVI.grad!(vo, alg, getq, model, param, diff_results[i], elbo_num)
        end

        ∇ = mean([DiffResults.gradient(diff_results[i]) for i in 1:minibatchsize])
        Δ = Flux.Optimise.apply!(optimizer, param, ∇)
        if any(isnan.(Δ))
            println("NaN occurred")
            continue
        end
        @. param = param - Δ

        ##########################################################
        push!(elbo_data, vo(alg, getq(param), model, elbo_num))

        if mod(i, 5) == 0.0 
            model   = fitLoss(X0[1], R[1], param; totalTimeStep=6000)
            converged = hasconverged(X0[1], R[1], param, elbo_data,i) 
        end

        ProgressMeter.next!(prog)
    end
end
