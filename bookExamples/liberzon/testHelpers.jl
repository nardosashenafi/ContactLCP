

function testBayesian(xi, par; totalTimeStep = totalTimeStep)
    θ   = rand(q(par))
    ψ, uk = unstackSamples(θ)

    testBayesian(xi, ψ, uk; totalTimeStep = totalTimeStep)
end

function testBayesian(xi, ψ, uk; totalTimeStep = totalTimeStep)    
    X,U = integrate(xi, ψ, uk; totalTimeStep = totalTimeStep)

    clf()
    plt.subplot(2, 2, 1)
    plot(getindex.(X, 1), getindex.(X, 2))
    scatter(X[end][1], X[end][2], marker="o")
    ylabel("x2")
    xlabel("x1")

    plt.subplot(2, 2, 2)
    plot(getindex.(X[1:end-1], 1), U)
    ylabel("uk")

    plotPartition(ψ, uk, X)
end

function testBayesianPaper(xi, ψ, uk; totalTimeStep = totalTimeStep)    
    X,U = integrate(xi, ψ, uk; totalTimeStep = totalTimeStep)

    clf()
    width = 30
    X1 = range(-10.0f0, 10.0f0, length=width)
    X2 = range(-10.0f0, 10.0f0, length=width)
    
    u = Matrix{Int}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            pk      = bin(vcat(X1[j], X2[width-i+1]), ψ)
            c[i, j] = argmax(pk)
            u[i, j] = uk[c[i, j]]
        end
    end

    custom_font = 15
    plt.subplot(1, 2, 1)
    scatter(X[end][1], X[end][2], marker="*", color="red", s=300, zorder=2)
    plot(getindex.(X, 1), getindex.(X, 2), color="black", linewidth=5, zorder=1)
    ylabel(L"x_2", fontsize=custom_font)
    xlabel(L"x_1", fontsize=custom_font)
    tick_params(axis="both", labelsize=custom_font)

    imshow(u, extent = [X1[1], X1[end], X2[1], X2[end]])
    title("Control", fontsize=custom_font)

    plt.subplot(1, 2, 2)
    imshow(c, extent = [X1[1], X1[end], X2[1], X2[end]])
    title("State partitions", fontsize=custom_font)
    ylabel(L"x_2", fontsize=custom_font)
    xlabel(L"x_1", fontsize=custom_font)
    tick_params(axis="both", labelsize=custom_font)
end

function plotPartition(ψ, uk, X)
    width = 30
    X1 = range(-5.0f0, 5.0f0, length=width)
    X2 = range(-5.0f0, 5.0f0, length=width)
    
    u = Matrix{Int}(undef, width, width)
    c = Matrix{Int}(undef, width, width)

    for i in 1:width    #row
        for j in 1:width    #column
            pk      = bin(vcat(X1[j], X2[width-i+1]), ψ)
            c[i, j] = argmax(pk)
            u[i, j] = uk[c[i, j]]
        end
    end

    plt.subplot(2, 2, 3)
    plot(getindex.(X, 1), getindex.(X, 2))
    scatter(X[end][1], X[end][2], marker="o")
    ylabel("x2")
    xlabel("x1")

    imshow(u, extent = [X1[1], X1[end], X2[1], X2[end]])
    ylabel("Control")

    plt.subplot(2, 2, 4)
    imshow(c, extent = [X1[1], X1[end], X2[1], X2[end]])
    ylabel("Partitions")
end


function checkGradient()

    mψ          = rand(Float32, nn_length)
    σψ          = LogExpFunctions.invsoftplus.(0.1 .* rand(Float32, nn_length))
    θk          = LogExpFunctions.logit.([0.8f0, 0.8f0, 0.3f0, 0.3f0])      
    vo          = Variational.ELBO()
    elbo_num    = 5
    alg         = ADVI(elbo_num, 100)

    par         = vcat(mψ, σψ, θk)
    diff_results = DiffResults.GradientResult(par)

    minibatch   = 2
    x0      = sampleInitialState(par; totalTimeStep=1000, minibatch=minibatch)

    model   = stateBinModel(x0, par; totalTimeStep=1000)
    Random.seed!(1234)
    AdvancedVI.grad!(vo, alg, q, model, par, diff_results, elbo_num)
    ∇       = DiffResults.gradient(diff_results)

    δ = 0.001
    par1 = deepcopy(par)
    par2 = deepcopy(par1)
    par2[1] += δ

    Random.seed!(1234)
    el1 = -vo(alg, q(par1), model, elbo_num)

    Random.seed!(1234)
    el2 = -vo(alg, q(par2), model, elbo_num)
    
    fd = abs(el2 - el1)/δ

    print("VI grad = ", ∇[1], " | fd = ", fd, " | error = ", abs(∇[1] - fd))
end

function elbo()

end

function evaluateTraining(chains, x0)
    b = [chains[["b[$i]"]].value.data[end] for i in 1:totalTimeStep]
    u = [chains[["u[$i]"]].value.data[end] for i in 1:binSize]
    t = [chains[["transition[$i]"]].value.data[end] for i in 1:binSize]
    θ = [chains[["θ[$i]"]].value.data[end] for i in 1:nn_length]

    X = trajectory(x0, p)

end

