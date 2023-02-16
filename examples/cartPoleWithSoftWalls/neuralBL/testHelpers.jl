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
    
    x1    = range(-2.0, 2.0, length=width)
    θdot     = range(-pi, pi , length=width)
    xdot = range(-0.5, 0.5, length=width)

    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(x1[j], 0.174, 0.09f0, θdot[width-i+1])
            pk      = bin(x, ψ)
            c[i,j]  = argmax(pk)
            u[i, j] = input(x, θk, c[i,j])[1]
        end
    end

    plt.subplot(2, 2, 2)
    imshow(c, extent = [x1[1], x1[end], θdot[1], θdot[end]])
    # contourf(x1, θ, u)
    ylabel("Partitions")
    xlabel("x1 vs θdot")
    # fig.colorbar(u, ax=ax1)
    #overlap the trajectory on the control heat map
    # plot(getindex.(X, 2), getindex.(X, 4), "r")
    # scatter(X[end][2], X[end][4])

    ################################################
    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(d+D, 0.174, xdot[j], θdot[width-i+1])
            pk      = bin(x, ψ)
            c[i,j]  = argmax(pk)
            u[i, j] = input(x, θk, c[i,j])[1]
        end
    end

    plt.subplot(2, 2, 4)
    imshow(c, extent = [xdot[1], xdot[end], θdot[1], θdot[end]])
    # myc = contourf(x1, θ, c, linestyles="None")
    # clabel(myc, fmt="%d")
    ylabel("State Partitions")
    xlabel("xdot vs θdot")

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

function activeGap(pendulum_xy, x1, θ1)

    if pendulum_xy[2] > rbl[2] || norm(pendulum_xy - rbl, 2) < gThreshold
        gi1,_, wi1,_  = gapPendulumToRightWall(pendulum_xy, θ1)
        gi2, wi2    = gapPendulumToBottomWall(pendulum_xy, rbl[2], θ1)
        gi3, gi4, wi3, wi4 = gapPendulumToRightWallCorners(pendulum_xy, x1, θ1)  

        gni = [gi1, gi2, gi3, gi4]
        _, ind = findmin(abs.(gni))
        gi = gni[ind]
        wi = [wi1, wi2, wi3, wi4][ind]
    else
        gi1 = norm(pendulum_xy - rbl, 2)
        x̄r, ȳr = pendulum_xy - rbl
        wi1 = 1.0f0/gi1*[x̄r; -x̄r*l*cos(θ1) + ȳr*l*sin(θ1)] 
        gi2, wi2 = gapPendulumToBottomWall(pendulum_xy, rbl[2], θ1)
        gni = [gi1, gi2]
        _, ind = findmin(abs.(gni))
        gi = gni[ind]
        wi = [wi1, wi2][ind]
    end
    return gi, wi
end

function gapvsRelativeVelvsControlPlot(ψ, θk)
    widthx1 = 15 
    widthθ1 = 15
    widthx1dot = 15 
    widthθ1dot = 15
    x1 = range(0.0, d+D, length=widthx1)
    θ1 = range(0.0, 2.0f0*pi, length=widthθ1)
    x1dot = range(-2.0f0, 2.0f0, length=widthx1dot)
    θ1dot = range(-14.0f0, 14.0f0, length=widthθ1dot)

    u_width = widthx1*widthθ1*widthx1dot*widthθ1dot
    u_height = widthx1*widthθ1*widthx1dot*widthθ1dot
    gn = Vector{Float32}(undef, u_width)
    gndot = Vector{Float32}(undef, u_height)
    u = Vector{Float32}(undef, u_width)
    c = Vector{Int}(undef, u_width)

    ind = 1

    for i in 1:widthx1    #row
        for j in 1:widthθ1    #column  
            for k in 1:widthx1dot    
                for m in 1:widthθ1dot   
                    x       = vcat(x1[i], θ1[j], x1dot[k], θ1dot[m])
                    pendulum_xy = pendulumPos(x1[i], θ1[j])
                    gi, wi = activeGap(pendulum_xy, x1[i], θ1[j])

                    gn[ind]    = gi
                    gndot[ind] = wi'*[x1dot[k], θ1dot[m]]
                    pk         = bin(x, ψ)
                    c[ind]     = argmax(pk)
                    u[ind]     = input(x, θk, c[ind])[1]
                    ind += 1
                end
            end
        end
    end

    return u, c, gn, gndot
end


function gapVsControlPlot(ψ, θk)
    widthx1 = 20 
    widthθ1 = 20
    widthx1dot = 1 
    widthθ1dot = 1
    x1 = range(0.0, d+D, length=widthx1)
    θ1 = range(0.0, 2*pi, length=widthθ1)
    x1dot = 0.0f0
    θ1dot = 0.0f0

    u_width = widthx1*widthθ1*widthx1dot*widthθ1dot
    gn = Matrix{Float32}(undef, widthx1, widthθ1)
    u = Matrix{Float32}(undef, widthx1, widthθ1)
    c = Matrix{Int}(undef, widthx1, widthθ1)

    for i in 1:widthx1    #row
        for j in 1:widthθ1    #column    
            x       = vcat(x1[j], θ1[i], x1dot, θ1dot)
            pendulum_xy = pendulumPos(x1[j], θ1[i])
            gi, wi = activeGap(pendulum_xy, x1[i], θ1[j])

            gn[i, j]    = gi
            pk          = bin(x, ψ)
            c[i, j]     = argmax(pk)
            u[i, j]     = input(x, θk, c[i, j])[1]
        end
    end

    return u, c, gn, x1, θ1
end

function discretizedStatesVsBin(ψ, θk)
    width =30
    θstar = -20.0f0*pi/180f0
    x1    = range(0.0, d+D, length=width)
    θdot  = range(-10.0f0, 10.0f0, length=width)
    # x1dot = l*cos(θstar) .* θdot
    x1dot = 0.0f0

    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(x1[j], θstar, x1dot, θdot[i])
            pk      = bin(x, ψ)
            c[i,j]  = argmax(pk)
            u[i, j] = input(x, θk, c[i,j])[1]
        end
    end

    return u, c, x1, θdot
end

function completePlots2(ψ, θk; xi = [0.0f0, pi, 0.0f0, 0.1f0])

    X, _ = integrate(xi, ψ, θk; totalTimeStep = 4500)
    fig1, ax1 = plt.subplots()
    # fig1 = figure()
    # ax1 = gca(projection="3d")
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    custom_fontsize=18
    custom_linewidth=2.0
    custom_colors = "black"
    #####################Plot Control input 

    #########x, θdot over the null space of Wn
    # u, c, x1, θdot = discretizedStatesVsBin(ψ, θk)
    # ax1.contourf(x1, θdot, c, cmap="PuBu")

    #######Scatter plot of g, ġ, bin
    # u, c, gn, gndot = gapvsRelativeVelvsControlPlot(ψ, θk)
    # ax1.scatter(gn, gndot, c=c)

    ########gap contour overlapped by bin lines
    u, c, gn, x1, θ1 = gapVsControlPlot(ψ, θk)
    gapContour = ax1.contourf(x1, θ1, gn, cmap="PuBu")
    binOngapContour = ax1.contour(x1, θ1, c, cmap="autumn")
    # clabel(binOngapContour, [1, 2], fmt="%d")
    cbar = colorbar(gapContour, ax = ax1)
    cbar.ax.tick_params(labelsize=custom_fontsize)
    ax1.set_title("Controller choice based on gap function at "* L"[\dot{x}, \dot{\theta}] = [0,0] ", fontsize=custom_fontsize)
    ax1.set_ylabel(L"\theta", fontsize=custom_fontsize)
    ax1.set_xlabel(L"x", fontsize=custom_fontsize)
    ax1.tick_params(axis="both", labelsize=custom_fontsize)

    #############post impact plots
    #################Plot trajectory
    ax2.plot(getindex.(X[1:3980], 2), getindex.(X[1:3980], 4), color=custom_colors, label="Simulation", linewidth=custom_linewidth)
    ax2.plot(getindex.(X[3980:4000], 2), getindex.(X[3980:4000], 4), linestyle="dotted", color=custom_colors, label="Impact", linewidth=custom_linewidth)
    ax2.plot(getindex.(X[4000:4500], 2), getindex.(X[4000:4500], 4), color=custom_colors)
    # ax2.scatter(X[3980][2], X[3980][4], marker="P", s=100, color="blue", zorder=3)
    ax2.arrow(X[750][2], X[750][4], 0.0, 0.5, head_width=0.1,head_length=1.0, fill=true, color=custom_colors)
    ax2.arrow(X[1300][2], X[1300][4], 0.0, -0.5, head_width=0.1,head_length=1.0, fill=true, color=custom_colors)
    # ax2.arrow(X[2100][2], X[2100][4], 0.0, 0.5, head_width=0.1,head_length=1.0, fill=true, color="red")
    # ax2.arrow(X[3100][2], X[3100][4], 0.0, -0.5, head_width=0.1,head_length=1.0, fill=true, color="red")
    ax2.arrow(X[3815][2], X[3815][4], -0.1, -0.0, head_width=0.5,head_length=0.15, fill=true, color=custom_colors)
    ax2.arrow(X[3982][2], X[3982][4]+3.0, 0.0, 0.5, head_width=0.1,head_length=1.0, fill=true, color=custom_colors)

    ##################plot experimental data
    file = MATLAB.MatFile("./hardware/cartpole_MOE/hardware_data/cartpole_MOE_2_8_firstCatchUnderlayment.mat", "r")
    data = get_variable(file, "cart_data")
    MATLAB.close(file)

    theta = 2*pi .- (pi .- data[4,:])
    thetadot = data[6,:]

    thetadot[1:100] .= 0.0
    ax2.plot(theta[1:3550], thetadot[1:3550],  linestyle="dashed", color=custom_colors, label="Experiment", linewidth=custom_linewidth)
    ax2.plot(theta[3550:3620], thetadot[3550:3620],  linestyle="dotted", color=custom_colors, linewidth=custom_linewidth)
    ax2.plot(theta[3620:3800], thetadot[3620:3800],  linestyle="dashed", color=custom_colors, linewidth=custom_linewidth)
    ax2.arrow(theta[650],  thetadot[650], -0.1, 0.0, head_width=0.5,head_length=0.15, fill=true, color=custom_colors)
    ax2.arrow(theta[1400], thetadot[1400], 0.0, -0.5, head_width=0.1,head_length=1.0, fill=true, color=custom_colors)
    # ax2.arrow(theta[2040], thetadot[2040], 0.0, 0.5, head_width=0.1,head_length=1.0, fill=true, color="green")
    # ax2.arrow(theta[3100], thetadot[3100], 0.0, -0.5, head_width=0.1,head_length=1.0, fill=true, color="green")
    ax2.arrow(theta[3450], thetadot[3450], -0.1, -0.0, head_width=0.5,head_length=0.15, fill=true, color=custom_colors)
    ax2.arrow(theta[3600]-0.01, thetadot[3600]-3.0, 0.0, 0.5, head_width=0.1,head_length=1.0, fill=true, color=custom_colors)
    ax2.scatter(theta[3800], thetadot[3800], marker="*", s=100, color="purple", zorder=2)

    ax2.legend(fontsize=custom_fontsize)

    #####################Plot Control input 

    X2    = range(-1.0, 2.0f0pi, length=width)
    X2dot = range(-14.0f0, 14.0f0, length=width)

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

    controlContour = ax2.contourf(X2, X2dot, u, cmap="PuBu")
    cbar2 = colorbar(controlContour, ax = ax2)
    cbar2.ax.tick_params(labelsize=custom_fontsize)
    ax2.set_title("Control Input at "* L"[x,\dot{x}]  = [0.36, 0.1]", fontsize=custom_fontsize)
    ax2.set_ylabel(L"\dot{\theta}", fontsize=custom_fontsize)
    ax2.set_xlabel(L"\theta", fontsize=custom_fontsize)
    ax2.tick_params(axis="both", labelsize=custom_fontsize)
    ###################Overlap bins contour plot
    binContour = ax2.contour(X2, X2dot, c, cmap="bone", zorder=1)
    # clabel(myc, fmt="%d")
    # ax2.scatter(X[3980][2], X[3980][4], marker="P", s=100, color="blue", zorder=2)
    # ax2.scatter(X[3990][2], X[3990][4], marker="*", s=100, color="red", zorder=2)
    # ax2.scatter(theta[3620], thetadot[3620], marker="*", s=100, color="green", zorder=2)

    #########show images

    postImpactImage = 1.0 .- Float64.(Gray.(load("../../media/plots/MOEpostImpact.png")))
    ax3.imshow(postImpactImage, cmap="binary")


end