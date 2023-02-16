using CSV 
using DataFrames 

function loadWeights(param)
    df =  DataFrame(param = param)
    CSV.write("neuralBL/savedWeights/hardwareParam2.csv", df)
end

function evaluateControl(x, ψ, θk)
    pk = bin(x, ψ)
    k = argmax(pk)
    u = input(x, θk, k)
    return [control(x, u)]
end 

function testBayesian(xi, par; totalTimeStep = totalTimeStep)
    θ    = rand(q(par))
    ψ, uk = unstackSamples(θ)

    testBayesian(xi, ψ, uk; totalTimeStep = totalTimeStep)
end

function testBayesian(xi, ψ, θk; totalTimeStep = totalTimeStep)    
    X, _ = integrate(xi, ψ, θk; totalTimeStep = totalTimeStep)
    clf()
    plots(X, fig1)
    # plotPartition2(X, ψ, θk)
    plotPartition(X, ψ, θk)
    fig1.canvas.draw()      #draws tupdates in for loop
    fig1.canvas.flush_events()  #gets new figure in for loop
    # animate(X)
    # sleep(1)
    return X
end

function plots(X)
    fig1 = plt.figure(1)
    plots(X, fig1)
end

function plots(X, fig1)
    plt.subplot(2, 2, 1)
    plot(getindex.(X, 2), getindex.(X, 4))
    scatter(X[end][2], X[end][4])
    ylabel("thetadot")
    xlabel("theta")

    plt.subplot(2, 2, 3)
    plot(getindex.(X, 1), getindex.(X, 3))
    scatter(X[end][1], X[end][3])
    ylabel("xdot")
    xlabel("x")
end


function plotPartition2(X::Vector{Vector{T}}, ψ, θk) where {T<:Real}
    width = 30
    
    X2    = range(-2.0f0pi, 2.0f0pi, length=width)
    X2dot = range(-10.0f0, 10.0f0, length=width)

    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(D+d, X2[j], 0.0f0, X2dot[i])
            pk      = bin(x, ψ)
            c[i,j]  = argmax(pk)
            u[i, j] = input(x, θk, c[i,j])[1]
        end
    end

    plt.subplot(2, 2, 2)
    # imshow(u, extent = [X2[1], X2[end], X2dot[1], X2dot[end]])
    contourf(X2, X2dot, u)
    ylabel("Control")
    xlabel("x2 vs x2dot")
    # fig.colorbar(u, ax=ax1)
    #overlap the trajectory on the control heat map
    # plot(getindex.(X, 2), getindex.(X, 4), "r")
    # scatter(X[end][2], X[end][4])

    ################################################
    plt.subplot(2, 2, 4)
    # imshow(c, extent = [X2[1], X2[end], X2dot[1], X2dot[end]])
    myc = contourf(X2, X2dot, c, linestyles="None")
    # clabel(myc, fmt="%d")
    ylabel("State Partitions")
    xlabel("x2 vs x2dot")

    ####################################################
    X1    = range(d+w/2, D+d-w/2, length=width)
    X1dot = range(d+w/2, D+d-w/2, length=width)

    u = Matrix{Float32}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(X1[j], 0.0f0, X1dot[i], 0.0f0)
            pk      = bin(x, ψ)
            c       = argmax(pk)
            u[i, j] = input(x, θk, c)[1]
        end
    end

    # plt.subplot(2, 2, 3)

    # imshow(u, extent = [X1[1], X1[end], X1dot[1], X1dot[end]])
    # ylabel("Control")
    # xlabel("x1 vs x1dot")

end

function plotPartition(X::Vector{Vector{T}}, ψ, θk) where {T<:Real}
    width = 30
    
    X2    = range(-2.0f0pi, 2.0f0pi, length=width)
    X2dot = range(-10.0f0, 10.0f0, length=width)

    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(d+D, X2[j], 0.0f0, X2dot[width-i+1])
            pk      = bin(x, ψ)
            c[i,j]  = argmax(pk)
            u[i, j] = input(x, θk, c[i,j])[1]
        end
    end

    plt.subplot(2, 2, 2)
    imshow(u, extent = [X2[1], X2[end], X2dot[1], X2dot[end]])
    ylabel("Control")
    xlabel("x2 vs x2dot")
    # fig.colorbar(u, ax=ax1)
    #overlap the trajectory on the control heat map
    plot(getindex.(X, 2), getindex.(X, 4), "r")
    scatter(X[end][2], X[end][4])

    ################################################
    plt.subplot(2, 2, 4)
    imshow(c, extent = [X2[1], X2[end], X2dot[1], X2dot[end]])
    ylabel("State Partitions")
    xlabel("x2 vs x2dot")

    ####################################################
    X1    = range(d+w/2, D+d-w/2, length=width)
    X1dot = range(d+w/2, D+d-w/2, length=width)

    u = Matrix{Float32}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(X1[j], 0.0f0, X1dot[width-i+1], 0.0f0)
            pk      = bin(x, ψ)
            c       = argmax(pk)
            u[i, j] = input(x, θk, c)[1]
        end
    end

    # plt.subplot(2, 2, 3)

    # imshow(u, extent = [X1[1], X1[end], X1dot[1], X1dot[end]])
    # ylabel("Control")
    # xlabel("x1 vs x1dot")

end

function plotPartition(ψ, θk, x1)
    width = 30
    
    X2    = range(-2.0f0pi, 2.0f0pi, length=width)
    X2dot = range(-10.0f0, 10.0f0, length=width)

    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(x1, X2[j], 0.0f0, X2dot[width-i+1])
            pk      = bin(x, ψ)
            c[i,j]  = argmax(pk)
            u[i, j] = input(x, θk, c[i,j])[1]
        end
    end

    plt.subplot(2, 2, 2)
    imshow(u, extent = [X2[1], X2[end], X2dot[1], X2dot[end]])
    ylabel("Control")
    xlabel("x2 vs x2dot")

    ################################################
    plt.subplot(2, 2, 4)
    imshow(c, extent = [X2[1], X2[end], X2dot[1], X2dot[end]])
    ylabel("State Partitions")
    xlabel("x2 vs x2dot")

end

function completePlots(ψ, θk; xi = [0.0f0, pi, 0.0f0, 0.1f0])

    X, _ = integrate(xi, ψ, θk; totalTimeStep = 4500)
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(figsize=(13, 3), ncols=3, nrows=2)
    ###################Plot trajectory
    
    ax1.plot(getindex.(X[1:3980], 2), getindex.(X[1:3980], 4), color="black", label="Simulation")
    ax1.scatter(X[1][2], X[1][4], marker="P", s=100, color="blue", zorder=2)
    ax1.scatter(X[3980][2], X[3980][4], marker="*", s=100, color="red", zorder=2)

    #################plot experimental data
    file = MATLAB.MatFile("./hardware/cartpole_MOE/hardware_data/cartpole_MOE_2_8_firstCatchUnderlayment.mat", "r")
    data = get_variable(file, "cart_data")
    MATLAB.close(file)

    theta = 2*pi .- (pi .- data[4,:])
    thetadot = data[6,:]

    thetadot[1:100] .= 0.0

    ax1.plot(theta[1:3550], thetadot[1:3550],  linestyle="dotted", color="black", label="Experiment")

    ax1.legend()
    ax1.set_ylabel(L"\dot{\theta}")
    ax1.set_xlabel(L"\theta")
    #####################Plot Control input 
    width = 30
    
    X2    = range(-1.0, 2.0f0pi, length=width)
    X2dot = range(-10.0f0, 10.0f0, length=width)

    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(0.0, X2[j], 0.0f0, X2dot[i])
            pk      = bin(x, ψ)
            c[i,j]  = argmax(pk)
            u[i, j] = input(x, θk, c[i,j])[1]
        end
    end

    controlContour = ax3.contourf(X2, X2dot, u, cmap="binary", zorder=1)
    colorbar(controlContour, ax = ax3)
    ax3.set_title("Control Input at "* L"[x,\dot{x}]  = [0, 0]")
    ax3.set_ylabel(L"\dot{\theta}")
    ax3.set_xlabel(L"\theta")

    ###################Overlap bins contour plot
    binContour = ax3.contour(X2, X2dot, c, cmap="bone")
    ax3.scatter(X[1][2], X[1][4], marker="P", s=100, color="blue", zorder=2)

    #############Repeat for post impact 
    #################Plot trajectory
    ax2.plot(getindex.(X[1:3980], 2), getindex.(X[1:3980], 4), color="black", label="Simulation")
    ax2.plot(getindex.(X[3980:4000], 2), getindex.(X[3980:4000], 4), linestyle="dashed", color="red", label="Impact")
    ax2.plot(getindex.(X[4000:4500], 2), getindex.(X[4000:4500], 4), color="black")
    ax2.scatter(X[3980][2], X[3980][4], marker="P", s=100, color="blue", zorder=2)
    ax2.scatter(X[3990][2], X[3990][4], marker="*", s=100, color="red", zorder=2)
    ax2.set_ylabel(L"\dot{\theta}")
    ax2.set_xlabel(L"\theta")

    ##################plot experimental data
    ax2.plot(theta[1:3550], thetadot[1:3550],  linestyle="dotted", color="black", label="Experiment")
    ax2.plot(theta[3550:3620], thetadot[3550:3620],  linestyle="dashed", color="green")
    ax2.scatter(theta[3620], thetadot[3620], marker="*", s=100, color="green", zorder=2)
    ax2.plot(theta[3620:3800], thetadot[3620:3800],  linestyle="dotted", color="black")
    ax2.legend()

    #####################Plot Control input 

    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(d+D, X2[j], 0.1f0, X2dot[i])
            pk      = bin(x, ψ)
            c[i,j]  = argmax(pk)
            u[i, j] = input(x, θk, c[i,j])[1]
        end
    end

    controlContour = ax4.contourf(X2, X2dot, u, cmap="binary")
    colorbar(controlContour, ax = ax4)
    ax4.set_title("Control Input at "* L"[x,\dot{x}]  = [0.36, 0.1]")
    ax4.set_ylabel(L"\dot{\theta}")
    ax4.set_xlabel(L"\theta")

    ###################Overlap bins contour plot
    binContour = ax4.contour(X2, X2dot, c, cmap="bone", zorder=1)
    # clabel(myc, fmt="%d")
    ax4.scatter(X[3980][2], X[3980][4], marker="P", s=100, color="blue", zorder=2)
    ax4.scatter(X[3990][2], X[3990][4], marker="*", s=100, color="red", zorder=2)
    ax4.scatter(theta[3620], thetadot[3620], marker="*", s=100, color="green", zorder=2)

    #########show images
    firstImpactImage = 1.0 .- Float64.(Gray.(load("../../media/plots/MOEfirstImpact.png")))
    ax5.imshow(firstImpactImage, cmap="binary")

    postImpactImage = 1.0 .- Float64.(Gray.(load("../../media/plots/MOEpostImpact.png")))
    ax6.imshow(postImpactImage, cmap="binary")


end