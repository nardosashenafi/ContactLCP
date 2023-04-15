include("EMSwingup.jl")

# oneController = BSON.load("./savedWeights/oneControllerHangingWalls_thetak5-8-8-4-4-1elu.bson")[:param]
moe = BSON.load("./savedWeights/for_paper_friction12_setdistancetraining_3NN_psi_5-4-4-3-3-1-thetak-5-10-10-4-4-1_elu copy.bson")[:param]
oneController = BSON.load("./savedWeights/oneController_friction12_setdistancetraining-thetak-5-10-10-4-4-1_elu.bson")[:param]
# nowalls = BSON.load("./savedWeights/swingupNoWallsThetak_5-8-8-4-4-1elu.bson")[:θk][1]

# function gap(x1, θ)
#     return Inf*ones(contactNum), Inf*ones(2, contactNum)
# end

# for i in 1:binSize 
#     controlArray[i] = FastChain(FastDense(5, 8, elu),
#                             FastDense(8, 4, elu),
#                             FastDense(4, 1))
    
#     controlNN_length[i] = DiffEqFlux.paramlength(controlArray[i]) 
# end

function oneTrajectory(xi, param::AbstractArray{T}; totalTimeStep = totalTimeStep) where {T<:Real}

    X = Vector{Vector{T}}(undef, totalTimeStep)
    X[1] = xi

    for i in 1:totalTimeStep-1
        u       = controlArray[1](inputLayer(X[i]), param)
        X[i+1]  = oneStep(X[i], u)
    end

    return X
end

function moeTraj(x0, param; totalTimeStep = totalTimeStep)
    ψ, θk   = unstackParams(param)
    X,_    = integrate(x0, ψ, θk; totalTimeStep = totalTimeStep)
    return X
end

function withinGoal(x)
    θ, θdot = [x[2], x[4]]
    cost = 1.0f0-cos(θ) + abs(θdot)

    return cost, ((1.0f0-cos(θ) <= (1.0f0-cosd(5.0))) && (abs(θdot) <= 0.3f0))
end 

function successCompare(oneController, moe; samples=100)

    moeSuccess = 0.0
    moeCost = 0.0
    oneControllerSuccess = 0.0
    oneControllerCost = 0.0

    x0MOEfail = []
    x0onefail = []
    for i in 1:samples
        x0  = Float32.([rand(-TRACK_LENGTH/2.0+w/2.0:0.01:TRACK_LENGTH/2.0-w/2.0), 
                            rand(0.0:0.005:2pi), 
                            rand(-1.0:0.001:1.0), 
                            rand(-3.0:0.001:3.0)])

        Xone = oneTrajectory(x0, oneController; totalTimeStep=10000);
        Xmoe = moeTraj(x0, moe; totalTimeStep=10000);

        onecost, onebool = withinGoal(Xone[end])
        moecost, moebool = withinGoal(Xmoe[end])

        if onebool
            oneControllerSuccess +=1
        else
            push!(x0onefail, x0)
        end
        oneControllerCost += onecost 
        if moebool
            moeSuccess +=1
        else
            push!(x0MOEfail, x0)
        end
        moeCost += moecost

    end
    return moeSuccess/samples, oneControllerSuccess/samples
end

function contourCompare(oneController, moe)
    x0 = [0.0f0, 3.1415f0, 0.0f0, -0.1f0]
    Xone = oneTrajectory(x0, oneController; totalTimeStep=9300);
    X1 = Xone[8850]
    Xmoe = moeTraj(Xone[7500], moe; totalTimeStep=2000);
    Xmoe1 = vcat(Xone[1:7500], Xmoe)


    clf()

    ulim = 10.0
    fig, (ax1, ax2) = plt.subplots(figsize=(15, 10), ncols=2, nrows=1)
    custom_linewidth= 4.0
    custom_fontsize = 25

    wid = 30
    X2    = range(-2.0, 2.0f0pi, length=wid)
    X2dot = range(-15.0f0, 13.0f0, length=wid)

    u = Matrix{Float32}(undef, wid, wid)

    for i in 1:wid    #row
        for j in 1:wid    #column
            x       = vcat(X1[1], X2[j], X1[3], X2dot[i])
            u[i, j] = clamp(controlArray[1](inputLayer(x), oneController)[1], -ulim, ulim)
        end
    end

    controlContour1 = ax1.contourf(X2, X2dot, u, cmap="PuBu", zorder=1, vmin=-ulim, vmax= ulim)
    # colorbar(controlContour1, ax = ax1)
    # ax1.set_title("Control Input at "* L"[x,\dot{x}]  = [0, 0]")
    ax1.set_ylabel(L"\dot{\theta}", fontsize=custom_fontsize)
    ax1.set_xlabel(L"\theta", fontsize=custom_fontsize)
    ax1.plot(getindex.(Xone[1:8850], 2), getindex.(Xone[1:8850], 4), color="black", linewidth=custom_linewidth)
    ax1.plot(getindex.(Xone[8850:8858], 2), getindex.(Xone[8850:8858], 4), linestyle="dashed", color="red", linewidth=custom_linewidth)
    ax1.plot(getindex.(Xone[8859:end], 2), getindex.(Xone[8859:end], 4), color="black", linewidth=custom_linewidth)
    ax1.scatter(Xone[end][2], Xone[end][4], marker="*", s=100, color="red", zorder=3)
    ax1.tick_params(axis="both", labelsize=custom_fontsize)

    ψ, θk = unstackParams(moe)
    um = Matrix{Float32}(undef, wid, wid)
    cm = Matrix{Int}(undef, wid, wid)

    for i in 1:wid    #row
        for j in 1:wid    #column
            x       = vcat(X1[1], X2[j], X1[3], X2dot[i])
            pk      = bin(x, ψ)
            cm[i,j]  = argmax(pk)
            um[i, j] = clamp(input(x, θk, cm[i,j])[1], -ulim, ulim)
        end
    end

    controlContour2 = ax2.contourf(X2, X2dot, um, cmap="PuBu", zorder=1, vmin=-ulim, vmax= ulim)
    cbar = colorbar(controlContour2, ax = ax2)
    cbar.ax.tick_params(labelsize=custom_fontsize)
    # ax2.set_title("Control Input at "* L"[x,\dot{x}]  = [0, 0]")
    ax2.set_ylabel(L"\dot{\theta}", fontsize=custom_fontsize)
    ax2.set_xlabel(L"\theta", fontsize=custom_fontsize)
    ax2.plot(getindex.(Xmoe1[1:7500], 2), getindex.(Xmoe1[1:7500], 4), color="black", linewidth=custom_linewidth, label="One controller")
    ax2.plot(getindex.(Xmoe1[7500:8630], 2), getindex.(Xmoe1[7500:8630], 4), linestyle="dashed", color="black", label="MOE", linewidth=custom_linewidth)
    ax2.plot(getindex.(Xmoe1[8630:8640], 2), getindex.(Xmoe1[8630:8640], 4), linestyle="dashed", color="red", label="Impact", linewidth=custom_linewidth)
    ax2.plot(getindex.(Xmoe1[8640:end], 2), getindex.(Xmoe1[8640:end], 4), color="black", linewidth=custom_linewidth)
    ax2.scatter(Xmoe1[end][2], Xmoe1[end][4], marker="*", s=100, color="red", zorder=3, label="Last state")

    binContour = ax2.contour(X2, X2dot, cm, cmap="bone")
    ax2.tick_params(axis="both", labelsize=custom_fontsize)
    handles, labels = ax2.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", fontsize=custom_fontsize, ncol=4)
    fig.tight_layout(pad=4.0)
end

function comparisonAnimation()
    x0 = [ 0.0
    3.1415
    0.067
    2.414]

    Xone = oneTrajectory(x0, nowalls; totalTimeStep=5000);
    Xmoe = moeTraj(x0, moe; totalTimeStep=7000);

end