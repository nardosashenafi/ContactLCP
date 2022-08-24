mutable struct Lcp{T, TSYS}
    sys                 ::TSYS 
    total_contact_num   ::Int
    current_contact_num ::Int
    ϵn                  ::Vector{T}
    ϵt                  ::Vector{T}
    μ                   ::Vector{T}
    gn                  ::Vector{T}
    γn                  ::Vector{T}
    γt                  ::Vector{T}
    ξn                  ::Vector{T}
    ξt                  ::Vector{T}
    Wn                  ::Matrix{T}
    Wt                  ::Matrix{T}

    function Lcp(sys, T)

        total_contact_num   = length(sys.contactIndex)
        current_contact_num = sum(sys.contactIndex)
        ϵn                  = sys.ϵn
        ϵt                  = sys.ϵt
        μ                   = sys.μ
        gn                  = zeros(T, total_contact_num)
        γn                  = zeros(T, total_contact_num)
        γt                  = zeros(T, total_contact_num)
        ξn                  = zeros(T, total_contact_num)
        ξt                  = zeros(T, total_contact_num)
        Wn                  = T.(wn(sys, sys.x0))
        Wt                  = T.(wt(sys, sys.x0)) 

        new{T, typeof(sys)}(sys, total_contact_num, current_contact_num, ϵn, ϵt, μ, gn, γn, γt, ξn, ξt, Wn, Wt)
    end
end

function checkContact(lcp::Lcp, gn)
     
    for i in 1:lcp.total_contact_num
        if gn[i] < G_THRESHOLD 
            lcp.sys.contactIndex[i] = 1 
        else
            lcp.sys.contactIndex[i] = 0 
        end
    end

    lcp.current_contact_num = sum(lcp.sys.contactIndex)
end

function resetConstants(lcp::Lcp)
    lcp.ϵn  = lcp.sys.ϵn
    lcp.ϵt  = lcp.sys.ϵt
    lcp.μ   = lcp.sys.μ
end

function createContactMap(lcp::Lcp, x, gn)

    lcp.gn = gn[lcp.sys.contactIndex .== 1]

    wn0 = wn(sys, x)
    wt0 = wt(sys, x)

    lcp.Wn = wn0[:, lcp.sys.contactIndex .== 1]   #pick out the ones in contact
    lcp.Wt = wt0[:, lcp.sys.contactIndex .== 1]

    lcp.γn = vnormal(lcp.sys, x)
    lcp.γt = vtang(lcp.sys, x)

    #trim the coefficients and relative velocities
    lcp.ϵn, lcp.ϵt, lcp.μ, lcp.γn, lcp.γt = map(x -> x[lcp.sys.contactIndex .== 1], 
                                            [lcp.ϵn, lcp.ϵt, lcp.μ, lcp.γn, lcp.γt])
end

function solveLcp(lcp::Lcp, M, h; Δt = 0.001)

    s   = lcp.current_contact_num

    if (s > 0)
        # println("Contact detected")
        E   = Matrix{Float64}(I, s, s)
        Wn  = lcp.Wn
        Wt  = lcp.Wt

        A = [Wn'*(M\(Wn - Wt*diagm(0 => lcp.μ))) Wn'*(M\Wt) zeros(s,s);
            Wt'*(M\(Wn - Wt*diagm(0 => lcp.μ))) Wt'*(M\Wt) E;
            2.0*diagm(0 => lcp.μ) -E zeros(s, s)]

        b = [Wn'*(M \ h*Δt) + (E + diagm(0 => lcp.ϵn))*lcp.γn;
            Wt'*(M \ h*Δt) + (E + diagm(0 => lcp.ϵt))*lcp.γt;
            zeros(s)]

        λ = lcpOpt(A, b, s)
    else
        λ = zeros(2s)
    end

    Λn = λ[1:s]
    ΛR = λ[s+1:2s]
    Λt = ΛR - diagm(0 => lcp.μ)*Λn
    
    return Λn, Λt, ΛR
end

function lcpOpt(A, b, contactNum)

    model = Model(PATHSolver.Optimizer)
    set_silent(model)
    # set_optimizer_attribute(model, "TimeLimit", 100)
    # set_optimizer_attribute(model, "OutputFlag", 1)
    # set_optimizer_attribute(model, "NonConvex", 2)
    @variable(model, λ[1:3*contactNum] >= 0.0)
    @constraints(model, begin
        (A*λ .+ b) ⟂ λ
    end)

    optimize!(model)
    # @show termination_status(model)


    return JuMP.value.(λ)
end


function onestep(lcp::Lcp, x1; Δt = 0.001)

    uA = x1[3:4]
    qA = x1[1:2]
    qM = qA + 0.5*Δt*uA

    x_mid = [qM...,uA...]
    resetConstants(lcp)
    gn = gap(lcp.sys, x_mid)

    checkContact(lcp, gn)
    createContactMap(lcp, x_mid, gn)

    M   = massMatrix(lcp.sys, x_mid)
    h   = genForces(lcp.sys, x_mid)

    Λn, Λt, ΛR = solveLcp(lcp, M, h; Δt=Δt)
    x2 = [qM...,uA...]

    λn = zeros(lcp.total_contact_num)
    λt = zeros(lcp.total_contact_num)
    λR = zeros(lcp.total_contact_num)
    λn[lcp.sys.contactIndex .== 1] = Λn
    λt[lcp.sys.contactIndex .== 1] = Λt
    λR[lcp.sys.contactIndex .== 1] = ΛR

    Wn = wn(lcp.sys, x_mid)
    Wt = wt(lcp.sys, x_mid)

    uE = M\((Wn - Wt*diagm(0 => lcp.sys.μ))*λn + Wt*λR + h*Δt) + uA
    qE = qM + 0.5*Δt*uE

    return [qE...,uE...], λn, λt

end

function fulltimestep(lcp::Lcp, T; Δt = 0.001)

    X       = Array{Array{T, 1}, 1}()
    Λn      = Array{Array{T, 1}, 1}()
    Λt      = Array{Array{T, 1}, 1}()
    gn      = Array{Array{T, 1}, 1}()
    t       = Array{T, 1}()
    x       = deepcopy(lcp.sys.x0)
    X       = push!(X, x)
    t       = push!(t, 0.0)

    for i in 1:1500
        x, λn, λt  = onestep(lcp, x; Δt=Δt)
        push!(X, x)
        push!(Λn, λn)
        push!(Λt, λt)
        push!(gn, gap(lcp.sys, x))
        push!(t, t[end]+Δt)
    end

    plots(X, t, Λn, Λt)
    return X, t, Λn, Λt, gn
end


function plots(Z, t, Λn, Λt)
    fig1 = plt.figure()
    fig1.clf()
    subplot(2, 2, 1)
    plot(t, getindex.(Z, 1))
    ylabel("x [m]", fontsize=15)
    subplot(2, 2, 2)
    plot(t, getindex.(Z, 2))
    ylabel("y [m]", fontsize=15)
    subplot(2, 2, 3)
    plot(t, getindex.(Z, 3))
    ylabel("vx [m/s]", fontsize=15)
    subplot(2, 2, 4)
    plot(t, getindex.(Z, 4))
    ylabel("vy [m/s]", fontsize=15)

    fig2 = plt.figure()
    fig2.clf()
    subplot(2, 1, 1)
    plot(t[2:end], getindex.(Λn, 1))
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{n1} [N]$", fontsize=15)
    # subplot(2, 3, 2)
    # plot(t[2:end], getindex.(Λn, 2))
    # ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    # ylabel(L"$\lambda_{n2} [N]$", fontsize=15)
    # subplot(2, 3, 3)
    # plot(t[2:end], getindex.(Λn, 3))
    # ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    # ylabel(L"$\lambda_{n3} [N]$", fontsize=15)
    subplot(2, 1, 2)
    plot(t[2:end], getindex.(Λt, 1))
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{t1} [N]$", fontsize=15)
    # subplot(2, 3, 5)
    # plot(t[2:end], getindex.(Λt, 2))
    # ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    # ylabel(L"$\lambda_{t2}$ [N]", fontsize=15)
    # subplot(2, 3, 6)
    # plot(t[2:end], getindex.(Λt, 3))
    # ticklabel_format(axis="y", style="sci",scilimits=(0,0))   
    # ylabel(L"$\lambda_{t3}$ [N]", fontsize=15)
end

# function animate(X)
    
#     @userplot ballPlot
#     @recipe function f(cp::ballPlot)
#         x, y, i = cp.args
#         n = length(x)
#         inds = circshift(1:n, 1 - i)
#         linewidth --> range(0, 10, length = n)
#         seriesalpha --> range(0, 1, length = n)
#         aspect_ratio --> 1
#         label --> false
#         x[inds], y[inds]
#     end

#     x = getindex.(X, 1)
#     y = getindex.(X, 2)
#     n = length(x)

#     anim = @animate for i ∈ 1:n
#         ballPlot(x, y, i)
#     end
#     gif(anim, "anim_fps15.gif", fps = 15)
# end