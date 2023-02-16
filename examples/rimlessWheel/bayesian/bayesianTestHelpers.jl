
include("bayesianTrainingHelpers.jl")

function testTrajectory(lcp, xi, sysParam, controlParam; Δt=0.0005f0, expert=false)
    X               = Vector{Vector{T}}(undef, totalTimeStep+1)
    Λn              = Vector{Vector{T}}(undef, totalTimeStep)
    ΛR              = Vector{Vector{T}}(undef, totalTimeStep)
    Λn_noise        = Vector{Vector{T}}(undef, totalTimeStep)
    ΛR_noise        = Vector{Vector{T}}(undef, totalTimeStep)
    X[1]            = deepcopy(xi)
 
    for i in 1:totalTimeStep
        X[i+1], Λn[i], ΛR[i], Λn_noise[i], ΛR_noise[i] = stateAndForcesWithNoise(lcp, X[i], sysParam, controlParam; Δt=Δt, expert=expert)
    end

    return X, Λn, ΛR, Λn_noise, ΛR_noise 
end

function checkContactForceNoise(lcp, xi, sysParam; controlParam=Float32[30.0, 5.0])
    X, Λn, ΛR, Λn_noise, ΛR_noise = testTrajectory(lcp, xi, sysParam, controlParam; Δt=0.0005f0, expert=true)

    t = range(0, (length(X)-1)*Δt, length=length(X)-1) 
    fig1 = figure()
    subplot(2, 2, 1)
    for i in 1:length(Λn[1])
        plot(t, getindex.(Λn, i))
    end
    ylabel(L"$\Lambda_n$")

    subplot(2, 2, 2)
    for i in 1:length(Λn_noise[1])
        plot(t, getindex.(Λn_noise, i))
    end
    ylabel(L"$\Lambda_n with noise$")

    subplot(2, 2, 3)
    for i in 1:length(ΛR[1])
        plot(t, getindex.(ΛR, i))
    end
    ylabel(L"$\Lambda_R$")

    subplot(2, 2, 4)
    for i in 1:length(ΛR_noise[1])
        plot(t, getindex.(ΛR_noise, i))
    end
    ylabel(L"$\Lambda_R with noise$")
end