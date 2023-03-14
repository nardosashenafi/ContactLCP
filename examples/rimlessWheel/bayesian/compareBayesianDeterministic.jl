responsetype    = DSP.Lowpass(0.6)
designmethod    = DSP.Butterworth(4)

function compareBayesianDeterministic(controlParam; k=k)

    Δh = range(0.0, step=0.1, length=10)

    for δh in Δh
        rvec = 0.0
        [append!(rvec[end] + iseven(i)*δh) for i in 1:k]
        initialStateWithBumps(θ0, θ0dot, ϕ0, ϕ0dot, rvec)
    end
end

function controllerDistributions(param, sampleNum)
    
    u   = Vector{eltype(param)}(undef, sampleNum)
    xi = initialState(3.1415f0, -0.5f0, 0.0f0, 0.0f0)

    for i in 1:sampleNum
        u[i] = MLBasedESC.controller(npbc, inputLayer(xi), rand(getq(param)))
    end

    Plots.histogram(u)
end

function wrapSpokes(θ)
    θnew  = θ
    if θ <= pi-α
        θnew = pi + α
    elseif θ >= pi + α
        θnew = pi - α
    end
    return θnew
end

function filter_signal(vel)
    DSP.filt(DSP.digitalfilter(responsetype, designmethod), vel)
end

function plotCollectedData(;δt=0.01)

    hardware_data = BSON.load("/home/nardosashenafi/repos/ContactLCP/examples/rimlessWheel/hardware_data/deterministic_bestRun.bson")[:sensorData]
    # hardware_data = BSON.load("/home/nardosashenafi/repos/ContactLCP/examples/rimlessWheel/hardware_data/bayesian_syncIssues.bson")[:sensorData]
    # hardware_data = BSON.load("/home/nardosashenafi/repos/ContactLCP/examples/rimlessWheel/hardware_data/b_gain_0_5.bson")[:sensorData]

    θdot    = getindex.(hardware_data, 4)
    ϕdot    = getindex.(hardware_data, 3)
    # θdot = getindex.(hardware_data, 6)
    θ       = getindex.(hardware_data, 2)
    xdot    = Vector{eltype(θdot)}(undef, length(θdot))
    t = range(1, stop=length(θdot)*δt, length=length(θdot))
    for i in eachindex(θdot)
        θIncontact  = wrapSpokes(θ[i])
        xdot[i]     = l1 * cos(θIncontact) * (θdot[i] + ϕdot[i])
    end
    
    # plot(t, filter_signal(xdot))
    plot(t, xdot)

    hardware_data = BSON.load("/home/nardosashenafi/repos/ContactLCP/examples/rimlessWheel/hardware_data/d_noGearRatio.bson")[:sensorData]
    θdot    = getindex.(hardware_data, 4)
    ϕdot    = getindex.(hardware_data, 3)
    # θdot2 = getindex.(hardware_data, 6)
    θ       = getindex.(hardware_data, 2)
    xdot    = Vector{eltype(θdot)}(undef, length(θdot))
    t = range(1, stop=length(θdot)*δt, length=length(θdot))

    for i in eachindex(θdot)
        θIncontact  = wrapSpokes(θ[i])
        xdot[i]     = l1 * cos(θIncontact) * (θdot[i]/6)
    end 
    
    # plot(t[end-10000:end], filter_signal(xdot[end-10000:end]))
    plot(t, xdot)
end