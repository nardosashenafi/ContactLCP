
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
const satu      = 2.0f0
responsetype    = DSP.Lowpass(0.5)
designmethod    = DSP.Butterworth(4)

function hipSpeedLoss(Z; gThreshold=gThreshold, k=k, α=α)

    #loss of one trajectory
    xd_dot  = 0.5f0

    loss = 0.0f0
    ki = 0

    for z in Z
        gn, _, _ = gap(z)
        contactIndex, _ = checkContact(z, gn, gThreshold, k)

        #loss between desired ẋ and actual ẋ should not be computed using getindex.(Z, 1) because this state does not directly depend on the control.
        #Using getindex.(Z, 1) as ẋ in auto-diff gives zero gradient
        if !isempty(contactIndex)
            ki = contactIndex[1] - 1
        else
            ki = spokeNearGround(gn) - 1
        end
        loss += xd_dot - (l1 * cos(z[4] + 2*α*ki) * z[8])

    end

    lmag = dot(loss, loss)

    ϕdot = getindex.(Z, 7)
    lmag += 2.0f0*dot(ϕdot, ϕdot)

    # #add cost on contact frequency
    # β       = 0.2f0
    # freq_d  = 1.0f0/(2.0f0*l1*sin(α)/xd_dot)

    # xdot          = getindex.(Z, 5)
    # θ             = getindex.(Z, 4)
    # xddot_avg     = 1.0f0/length(θ)*sum(xdot)
    # strikePeriod  = 2.0f0*l1*sin(α)/xddot_avg  #TODO: assumes  if the spoke in contact
    # f             = 1.0f0/strikePeriod
    # if f >= (1+β)*freq_d
    #     lmag += 5.0f0*(f - (1+β) .* freq_d)
    # end

    return 1.0f0/length(Z)*lmag
end

function extractStumbling(X)
    ind = findfirst(x -> x >= -0.01, getindex.(X, 8))
    if !isnothing(ind)
        return X[1:ind] 
    else
        return X
    end
end

@model function fitLoss(X0, R, μ_prior::AbstractArray{T}, σ_prior::AbstractArray{T}; totalTimeStep = 2000) where {T<:Real}

    # w  = Vector{T}(undef, length(μ_prior))
    # w .~ Distributions.Normal.(μ_prior, LogExpFunctions.softplus.(σ_prior))
    w ~ arraydist(map((m, s) -> Distributions.Normal(m, LogExpFunctions.softplus(s)), μ_prior, σ_prior))

    l11 = 0.0f0
    for i in eachindex(X0)
        X ,_  = trajectory(X0[i], R[i], w; totalTimeStep=totalTimeStep)  
        # X1  = extractStumbling(X)
        l11 += hipSpeedLoss(X)
    end
    
    # 0.0 ~ Normal(l11/batch_size, LogExpFunctions.softplus.(s[1]))
    0.0 ~ Distributions.Normal(l11/length(X0), 0.2f0)
end

function getq(param)
    μ_param, σ_param = unstackParams(param)
    return MvNormal(μ_param, LogExpFunctions.softplus.(σ_param))
    # return arraydist(map((m, σ) -> Distributions.Normal(m, LogExpFunctions.softplus.(σ)), μ_param, σ_param))
end

function unstackParams(param)
    μ_param = @view param[1:paramNum]
    σ_param = @view param[paramNum+1:end]

    return μ_param, σ_param
end

function controlToHipSpeed()

    fig1 = figure()
    fig2, ax2 = subplots()
    Turing.setrdcache(true)
    Turing.setadbackend(:forwarddiff)

    num_steps   = 5000
    elbo_num    = 5
    minibatchsize  = 2

    param = Float32.(vcat(0.3f0*randn(paramNum), invsoftplus.(0.1f0*rand(paramNum))))
    μ_param, σ_param = unstackParams(param)

    optimizer   = Variational.DecayedADAGrad(0.001)

    vo          = Variational.ELBO()
    alg         = ADVI(elbo_num, num_steps)
    converged   = false

    prog        = ProgressMeter.Progress(num_steps, 1)
    diff_results = DiffResults.GradientResult(param)
    elbo_data   = Vector{Float32}()

    @showprogress for i in 1:5000

        μ_param, σ_param = unstackParams(param)
        X0, R      = sampleInitialStates(param, minibatchsize; totalTime=3000)
        model   = fitLoss(X0, R, μ_param, σ_param; totalTimeStep=2000)
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

        ##########################################################
        mod(i, 10) == 0.0 ? converged = hasconverged(X0[1], param, elbo_data,i) : nothing

        ProgressMeter.next!(prog)
    end
end

function hasconverged(x0, param, elbo_data, i)

    filtered_elbo = filter_elbo(elbo_data[end-9:end])

    println("elbo = ", round(elbo_data[end], digits=4))

    w                   = rand(getq(param))
    X, sysParam         = trajectory(x0, [], w; totalTimeStep=10000)
    loss                = hipSpeedLoss(X)
    plots(X, fig1)
    ax2.plot(i, filtered_elbo[end], marker=".", color="k") 
    BSON.@save "./saved_weights/RW_bayesian_6-8-8-5-5-1_elu.bson" param
    println("loss = ", round(loss, digits=4) , " | hip speed = ", round.(mean(getindex.(X, 5)), digits=4) )

    return false
end

function filter_elbo(Δelbo_vec)
    DSP.filt(DSP.digitalfilter(responsetype, designmethod), Δelbo_vec)
end