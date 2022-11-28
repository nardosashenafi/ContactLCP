

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

    plt.subplot(2, 2, 2)
    plot(-l .* sin.(getindex.(X, 2)), l .* cos.(getindex.(X, 2)))
    scatter(l .* sin.(X[end][2]), l .* cos.(X[end][2]))
    ylabel("cos(theta)")
    xlabel("sin(theta)")
end

function plotPartition(ψ, θk)
    width = 30
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

