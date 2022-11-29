

function testBayesian(xi, par; totalTimeStep = totalTimeStep)
    θ    = rand(q(par))
    ψ, uk = unstackSamples(θ)

    testBayesian(xi, ψ, uk; totalTimeStep = totalTimeStep)
end

function testBayesian(xi, ψ, θk; totalTimeStep = totalTimeStep)    
    X, U = integrate(xi, ψ, θk; totalTimeStep = totalTimeStep)
    clf()
    plots(X, fig1)
    plotPartition(ψ, θk)
    animate(X)
    sleep(1)
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

end

function plotPartition(ψ, θk)
    width = 30
    
    X2    = range(-pi, pi, length=width)
    X2dot = range(-5.0, 5.0, length=width)

    u = Matrix{Float32}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(0.0f0, X2[j], 0.0f0, X2dot[width-i+1])
            pk      = bin(x, ψ)
            c       = argmax(pk)
            u[i, j] = input(x, θk, c)[1]
        end
    end

    plt.subplot(2, 2, 2)
    imshow(u, extent = [X2[1], X2[end], X2dot[1], X2dot[end]])
    ylabel("Control")
    xlabel("x2 vs x2dot")

    X1    = range(d, -d, length=width)
    X1dot = range(d, -d, length=width)
    
    u = Matrix{Float32}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            x       = vcat(X1[j], 0.0f0, X1dot[width-i+1], 0.0f0)
            pk      = bin(x, ψ)
            c[i, j] = argmax(pk)
            u[i, j] = input(x, θk, c[i, j])[1]
        end
    end


    plt.subplot(2, 2, 3)

    imshow(u, extent = [X1[1], X1[end], X1dot[1], X1dot[end]])
    ylabel("Control")
    xlabel("x1 vs x1dot")

    plt.subplot(2, 2, 4)
    imshow(c, extent = [X1[1], X1[end], X1dot[1], X1dot[end]])
    ylabel("Partitions")
    xlabel("x1 vs x1dot")

end

