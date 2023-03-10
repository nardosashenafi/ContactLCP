
function compareBayesianDeterministic(controlParam; k=k)

    Δh = range(0.0, step=0.1, length=10)
    hardware_data = BSON.load("/home/nardosashenafi/repos/ContactLCP/examples/rimlessWheel/hardware_data/determinisitic_bestRun.bson")

    for δh in Δh
        rvec = 0.0
        [append!(rvec[end] + iseven(i)*δh) for i in 1:k]
        initialStateWithBumps(θ0, θ0dot, ϕ0, ϕ0dot, rvec)
    end
end

function controllerDistributions(param, sampleNum)
    
    u   = Vector{eltype(param)}(undef, sampleNum)
    xi  = [2.782345f-5
             0.32197183
             0.0
             3.1415
             0.15
            -1.3911725f-5
             0.0
            -0.5]

    for i in 1:sampleNum
        u[i] = MLBasedESC.controller(npbc, inputLayer(xi), rand(getq(param)))
    end

    Plots.histogram(u)
end