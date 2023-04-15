include("bayesianControl.jl")
using CSV

Hd              = FastChain(FastDense(6, 8, elu), 
                  FastDense(8, 7, elu),
                  FastDense(7, 1))

npbc            = MLBasedESC.NeuralPBC(N, Hd)
const paramNum  = DiffEqFlux.paramlength(Hd)+N

responsetype    = DSP.Lowpass(0.6)
designmethod    = DSP.Butterworth(4)


function testLoss(Z, obstacles; gThreshold=gThreshold, k=k, α=α)

    #loss of one trajectory
    xd_dot  = 1.0f0
    loss = 0.0f0
    ki = 0

    for i in eachindex(Z)
        z = Z[i]
        gn, _, _ = gap(z, obstacles[i])
        contactIndex, _ = checkContact(z, gn, gThreshold, k)

        #loss between desired ẋ and actual ẋ should not be computed using getindex.(Z, 1) because this state does not directly depend on the control.
        #Using getindex.(Z, 1) as ẋ in auto-diff gives zero gradient
        if !isempty(contactIndex)
            ki = contactIndex[1] - 1
        else
            ki = spokeNearGround(gn) - 1
        end
        error = xd_dot - (l1 * cos(z[4] + 2*α*ki) * z[8])
        loss += 50.0f0*dot(error, error) #+ 0.5f0*dot(z[7], z[7])
    end

    return 1.0f0/length(Z)*loss
end

function deterTrajectory(x0, r, controlParam::Vector{T}; expert=false, Δt = Δt, totalTimeStep = 1000) where {T<:Real}

    X           = Vector{Vector{T}}(undef, totalTimeStep)
    obstacles   = Vector{Vector{Vector{T}}}(undef, totalTimeStep)
    x           = deepcopy(x0)
    sysParam    = createUnevenTerrain(x, r)
    θi          = x[4]

    for i in 1:totalTimeStep

        if abs(x[4] - θi) > 2pi - 2α
            sysParam = createUnevenTerrain(x, r)
            θi = x[4]
        end
        x       = oneStep(x, sysParam, controlParam; Δt=Δt, expert=expert)
        X[i]    = x
        obstacles[i] = sysParam
    end

    return X, obstacles
end

bayesMarginalTrajectory(x0, r, controlParam, sampleNum;expert=false, Δt = Δt, totalTimeStep = 1000) = integrateMarginalization(x0, r, controlParam, sampleNum; expert=expert, Δt = Δt, totalTimeStep = totalTimeStep) 

function bayesMapTrajectory(x0, r, controlParam::Vector{T}; expert=false, Δt = Δt, totalTimeStep = 1000) where {T<:Real}

    X           = Vector{Vector{T}}(undef, totalTimeStep)
    obstacles   = Vector{Vector{Vector{T}}}(undef, totalTimeStep)
    x           = deepcopy(x0)
    sysParam    = createUnevenTerrain(x, r)
    θi          = x[4]

    for i in 1:totalTimeStep

        if abs(x[4] - θi) > 2pi - 2α
            sysParam = createUnevenTerrain(x, r)
            θi = x[4]
        end
        u       = map(x, controlParam)
        x       = oneStep(x, sysParam, controlParam; Δt=Δt, expert=expert)
        X[i]    = x
        obstacles[i] = sysParam
    end

    return X, obstacles
end

###both trainings use the same size neural network
deterParam = BSON.load("../deterministic/saved_weights/deter2_hardware_even_1mpers_6-8-8-7-7-1_elu.bson")[:param]
# bayesParam = BSON.load("./saved_weights/RW_bayesian_6-8-8-5-5-1_elu.bson")[:param]
bayesParam = BSON.load("./saved_weights/RW_bayesian_6-8-8-7-7-1_elu.bson")[:param]

function compareBayesianDeterministic(deterParam, bayesParam; samples=10, k=k)
    rmax = [0.0f0, 0.005f0, 0.01f0, 0.015f0, 0.02f0]
    loss = []

    for rm in rmax
        l_deter = []
        l_marginal = []

        for i in 1:samples
            xi, r = initialStateWithBumps(
                        3.1415f0, 
                        -0.5f0, 
                        0.0f0, 
                        0.0f0, rm)
            X_deter, obstacle_deter = deterTrajectory(xi, r, deterParam; totalTimeStep = 20000);
            push!(l_deter, testLoss(X_deter, obstacle_deter))

            X_marginal, obstacle_marginal = bayesMarginalTrajectory(xi, r, bayesParam, 15; totalTimeStep = 20000);
            push!(l_marginal, testLoss(X_marginal, obstacle_marginal))
            println(i, "th iteration")
            # X_map, _ = bayesMapTrajectory(xi, r, bayesParam; totalTimeStep = 15000)
        end
        push!(loss, [l_deter, l_marginal])
    end
    return loss
end

function plotHistogramComparison(lossdata)

    loss0_deter         = lossdata[1][1]
    loss0_bayes         = lossdata[1][2]
    loss0_5_deter       = lossdata[2][1]
    loss0_5_bayes       = lossdata[2][2]
    loss1_deter         = lossdata[3][1]
    loss1_bayes         = lossdata[3][2]
    loss1_5_deter       = lossdata[4][1]
    loss1_5_bayes       = lossdata[4][2]
    loss2_deter         = lossdata[5][1]
    loss2_bayes         = lossdata[5][2]

    x = [0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.5, 1.5, 2.0, 2.0]
    height = [mean(loss0_deter), mean(loss0_bayes),
                mean(loss0_5_deter), mean(loss0_5_bayes),
                mean(loss1_deter), mean(loss1_bayes),
                mean(loss1_5_deter), mean(loss1_5_bayes),
                mean(loss2_deter), mean(loss2_bayes)]
    std_h = [std(loss0_deter), std(loss0_bayes),
                std(loss0_5_deter),std(loss0_5_bayes),
                std(loss1_deter), std(loss1_bayes),
                std(loss1_5_deter),std(loss1_5_bayes),
                std(loss2_deter), std(loss2_bayes)]

    x_deter = x[1:2:end]
    height_deter = height[1:2:end]
    x_bayes = x[2:2:end]
    height_bayes  = height[2:2:end]

    custom_fontsize = 25
    PyPlot.bar(x_deter, height_deter, width=0.3, color="blue", label="Deterministic")
    PyPlot.bar(x_bayes, height_bayes, width=0.3, color="orange", label="Bayesian")
    legend(fontsize=23)
    tick_params(axis="both", labelsize=custom_fontsize)
    xlabel(L"p_{max}", fontsize=custom_fontsize)
    ylabel(L"J_T", fontsize=custom_fontsize)
    ylim(top=1.3)
end

function controllerDistributions(param, sampleNum)
    
    u   = Vector{eltype(param)}(undef, sampleNum)
    xi = initialState(3.1415f0, -0.5f0, 0.0f0, 0.0f0)

    for i in 1:sampleNum
        u[i] = MLBasedESC.controller(npbc, inputLayer(xi), rand(getq(param)))
    end

    Plots.histogram(u)
end


function deterministicPlot(deterParam)
    custom_fontsize=18
    xi, r = initialStateWithBumps(
                        3.1415f0, 
                        -0.5f0, 
                        0.0f0, 
                        0.0f0, 0.0f0)
    Z, _ = deterTrajectory(xi, r, deterParam; totalTimeStep = 50000);
    t = range(0, step=Δt, length=length(Z))
    PyPlot.figure(1)
    fig1.clf()
    subplot(2, 1, 1)
    ϕdot = []
    ϕ = []
    push!(ϕdot, [Z[1][7]])
    push!(ϕ, [Z[1][3]])

    for i in 2:length(Z) 
        if abs(Z[i][7] - Z[i-1][7]) > 1.0
            push!(ϕdot,  [Z[i][7]])
            push!(ϕ,  [Z[i][3]])
        else
            push!(ϕdot[end], Z[i][7])
            push!(ϕ[end], Z[i][3])
        end
    end

    for i in eachindex(ϕdot)
        plot(ϕ[i], ϕdot[i], color="black")
    end
    for i in 1:length(ϕdot)-1
        plot([ϕ[i][end], ϕ[i+1][1]], [ϕdot[i][end], ϕdot[i+1][1]], color="red", "--")
    end

    ylabel(L"\dot{\phi} [rad/s]", fontsize=custom_fontsize)
    xlabel(L"{\phi} [rad]", fontsize=custom_fontsize)
    arrow(ϕ[1][265], ϕdot[1][265], 0.1, -0.0, width= 0.001, head_width=0.6,head_length=0.02, fill=true, color="black", length_includes_head=true)
    tick_params(axis="both", labelsize=custom_fontsize)
    arrow(1.0, 0.0, 0.0, -2.0, head_width=0.02,head_length=1.0, fill=true, color="red",  length_includes_head=true)

    subplot(2, 1, 2)
    plot(t, getindex.(Z, 5), color="black")
    ylabel(L"\dot{x} [m/s]", fontsize=custom_fontsize)
    xlabel(L"t [seconds]", fontsize=custom_fontsize)
    tick_params(axis="both", labelsize=custom_fontsize)

    ###############contour plot 
    
    width=50
    custom_fontsize=35
    θdot = range(-3.0, stop=2.5, length=width)
    xdot = -l1 .* θdot
    ϕ = range(0.0, stop=pi/2, length=width)
    u = Matrix{Float32}(undef, width, width)

    for i in eachindex(ϕ)
        for j in eachindex(θdot)
            x = [cos(ϕ[i]), sin(ϕ[i]), cos(pi), sin(pi), 0.0f0, θdot[j]]
            u[j, i] = clamp(MLBasedESC.controller(npbc, x, deterParam), -satu, satu)
        end
    end

    fig, (ax1) = plt.subplots(figsize=(13, 3), ncols=1, nrows=1)

    controlContour = ax1.contourf(ϕ, xdot, u, cmap="binary", zorder=1)
    cbar = colorbar(controlContour, ax = ax1)
    cbar.ax.tick_params(labelsize=custom_fontsize)
    ax1.set_title("Control Input at "* L"[\theta,\dot{\phi}]  = [\pi, 0]",fontsize=custom_fontsize)
    ax1.set_ylabel(L"\dot{x}", fontsize=custom_fontsize)
    ax1.set_xlabel(L"\phi", fontsize=custom_fontsize)
    ax1.tick_params(axis="both", labelsize=custom_fontsize)
end


function wrapSpokes!(θnew, i)
    if θnew[i] <= pi-α
        θnew[i:end] .+= pi + α
    elseif θnew[i] >= pi + α
        θnew[i:end] .+= pi - α
    end
end

function wrapSpokes(θ)
    θnew = θ
    if θnew <= pi-α
        θnew = pi + α
    elseif θnew >= pi + α
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
        xdot[i]     = l1 * cos(θIncontact) * (θdot[i])
    end
    
    # plot(t, filter_signal(xdot))
    plot(t, xdot)

    # hardware_data = BSON.load("/home/nardosashenafi/repos/ContactLCP/examples/rimlessWheel/hardware_data/d_even.bson")[:sensorData]
    hardware_data = BSON.load("/home/nardosashenafi/repos/ContactLCP/examples/rimlessWheel/hardware_data/d_nogainNogearRation.bson")[:sensorData]

    θdot    = getindex.(hardware_data, 4) 
    ϕdot    = getindex.(hardware_data, 3)
    ϕ       = getindex.(hardware_data, 1)
    # θdot = getindex.(hardware_data, 6)
    θ       = getindex.(hardware_data, 2) 
    xdot    = Vector{eltype(θdot)}(undef, length(θdot))
    t       = range(1, stop=length(θdot)*δt, length=length(θdot))
    θInContact = deepcopy(θ)

    for i in eachindex(θdot)
        θInContact = wrapSpokes(θ[i])
        xdot[i]     = l1 * cos(θInContact) * (θdot[i]/4)
    end
    
    plot(t, xdot)
    plot(t[66150:66850], xdot[66150:66850])
end

function plotRosbagDeterUneven()

    filethetadot = "./../hardware_data/rosbags/thetadot_gravel1.csv"
    filetheta = "./../hardware_data/rosbags/theta_gravel1.csv"

    csv_readerthetadot=CSV.File(filethetadot; delim="---")
    csv_readertheta=CSV.File(filetheta; delim="---")
    thetadot = []
    theta = []

    for i in 2:2:length(csv_readerthetadot)
        push!(thetadot, parse(Float32, csv_readerthetadot[i][1]))
    end
    for i in 2:2:length(csv_readertheta)
        push!(theta, parse(Float32, csv_readertheta[i][1]))
    end
    xdot = []
    for i in eachindex(thetadot)
        θInContact = wrapSpokes(pi + theta[i])
        push!(xdot, l1 * cos(θInContact) * (thetadot[i])*1/0.8)
    end

    custom_linewidth= 2.0
    custom_fontsize = 25
    t = range(1, length=length(thetadot), step=1/40)
    start = 1
    stop=800
    plot(range(0, stop=t[stop], length=length(t[start:stop])), xdot[start:stop], linewidth=custom_linewidth, color="black")
    plot(range(0, stop=t[stop], length=length(t[start:stop])), ones(length(t[start:stop])), color="red", linewidth=2*custom_linewidth, linestyle="dashed", label="Desired Hip Speed")


    ylabel(L"\dot{x} [m/s]", fontsize = custom_fontsize)
    xlabel(L"t [seconds]", fontsize = custom_fontsize)
    tick_params(axis="both", labelsize=custom_fontsize)
    legend(fontsize = custom_fontsize)
end

function plotRosbagDeterEven()

    filethetadot = "./../hardware_data/rosbags/thetadot_labLevel2.csv"
    filetheta = "./../hardware_data/rosbags/theta_labLevel2.csv"

    # filethetadot = "./../hardware_data/rosbags/thetadot_labLevel3.csv"
    # filetheta = "./../hardware_data/rosbags/theta_labLevel3.csv"

    csv_readerthetadot=CSV.File(filethetadot; delim="---")
    csv_readertheta=CSV.File(filetheta; delim="---")
    thetadot = []
    theta = []

    for i in 2:2:length(csv_readerthetadot)
        push!(thetadot, parse(Float32, csv_readerthetadot[i][1]))
    end
    for i in 2:2:length(csv_readertheta)
        push!(theta, parse(Float32, csv_readertheta[i][1]))
    end
    xdot = []
    for i in eachindex(thetadot)
        θInContact = wrapSpokes(pi + theta[i])
        push!(xdot, l1 * cos(θInContact) * (thetadot[i])*1/0.8)
    end

    custom_linewidth= 2.0
    custom_fontsize = 25
    t = range(1, length=length(thetadot), step=1/40)
    # plot(t, xdot, color="black", label="Deterministic on Even Terrain")
    start = 200
    stop=650
    plot(range(0, stop=t[stop], length=length(t[start:stop])), xdot[start:stop], linewidth=custom_linewidth, color="black")
    plot(range(0, stop=t[stop], length=length(t[start:stop])), ones(length(t[start:stop])), color="red", linewidth=2*custom_linewidth, linestyle="dashed", label="Desired Hip Speed")


    ylabel(L"\dot{x} [m/s]", fontsize = custom_fontsize)
    xlabel(L"t [seconds]", fontsize = custom_fontsize)
    tick_params(axis="both", labelsize=custom_fontsize)
    legend(fontsize = custom_fontsize)
end


function plotRosbagBayesEven()

    filethetadot = "./../hardware_data/rosbags/thetadot_bayes_levelLab1.csv"
    filetheta = "./../hardware_data/rosbags/theta_bayes_levelLab1.csv"

    # filethetadot = "./../hardware_data/rosbags/thetadot_labLevel3.csv"
    # filetheta = "./../hardware_data/rosbags/theta_labLevel3.csv"

    csv_readerthetadot=CSV.File(filethetadot; delim="---")
    csv_readertheta=CSV.File(filetheta; delim="---")
    thetadot = []
    theta = []

    for i in 2:2:length(csv_readerthetadot)
        push!(thetadot, parse(Float32, csv_readerthetadot[i][1]))
    end
    for i in 2:2:length(csv_readertheta)
        push!(theta, parse(Float32, csv_readertheta[i][1]))
    end
    xdot = []
    for i in eachindex(thetadot)
        θInContact = wrapSpokes(pi + theta[i])
        push!(xdot, l1 * cos(θInContact) * (thetadot[i])*1/0.8)
    end

    custom_linewidth= 2.0
    custom_fontsize = 25
    t = range(1, length=length(thetadot), step=1/40)
    # plot(t, xdot, color="black", label="Bayesian on Even Terrain")
    # plot(t, xdot, color="black")
    # plot(t, ones(length(t)), color="red", linewidth=custom_linewidth, linestyle="dashed", label="Desired Hip Speed")
    start = 1300
    stop=1800
    plot(range(0, stop=t[stop], length=length(t[start:stop])), xdot[start:stop], linewidth=custom_linewidth, color="black")
    plot(range(0, stop=t[stop], length=length(t[start:stop])), ones(length(t[start:stop])), color="red", linewidth=2*custom_linewidth, linestyle="dashed", label="Desired Hip Speed")


    ylabel(L"\dot{x} [m/s]", fontsize = custom_fontsize)
    xlabel(L"t [seconds]", fontsize = custom_fontsize)
    tick_params(axis="both", labelsize=custom_fontsize)
    legend(fontsize = custom_fontsize)
end