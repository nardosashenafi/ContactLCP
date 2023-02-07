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

function plotPartition(X::Vector{Vector{T}}, ψ, θk) where {T<:Real}
    width = 30
    
    X2    = range(-2.0f0pi, 2.0f0pi, length=width)
    X2dot = range(-10.0f0, 10.0f0, length=width)

    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(D+d, X2[j], 0.0f0, X2dot[width-i+1])
            pk      = bin(x, ψ)
            c[i,j]  = argmax(pk)
            u[i, j] = input(x, θk, c[i,j])[1]
        end
    end

    plt.subplot(2, 2, 2)
    imshow(u, extent = [X2[1], X2[end], X2dot[1], X2dot[end]])
    ylabel("Control")
    xlabel("x2 vs x2dot")

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