export BouncingBall

struct BouncingBall{T} 
    m               ::T 
    r               ::T
    g               ::T 
    ϵn              ::Vector{T}
    ϵt              ::Vector{T}
    μ               ::Vector{T}
    x0              ::Vector{T}
    contactIndex    ::Vector{T}
    gThreshold      ::T

    function BouncingBall(T)

        m               = T(0.2)
        r               = T(0.1)
        g               = T(9.81)
        ϵn              = T.(0.9*ones(1))
        ϵt              = T.(-0.5*ones(1))
        μ               = T.(0.2*ones(1))
        x0              = T.([0.1, 0.3, 0.1, -0.1])     #xpos, ypos, xdot, ydot
        contactIndex    = zeros(T, 1)
        gThreshold      = T.(0.001)

        new{T}(m, r, g, ϵn, ϵt, μ, x0, contactIndex, gThreshold)
    end

end

function (sys::BouncingBall)(x, param)
    gn  = gap(sys, x)
    γn  = vnormal(sys, x)
    γt  = vtang(sys, x)
    M   = massMatrix(sys, x)
    h   = genForces(sys, x)
    Wn  = wn(sys, x)
    Wt  = wt(sys, x)

    return gn, γn, γt, M, h, Wn, Wt
end

function (sys::BouncingBall)(x::Vector{T}) where {T<:Real}
    q = x[1:2]
    u = x[3:4]
    return q, u
end

function gap(sys::BouncingBall, x)
    return [x[2] - sys.r]
end

function vnormal(sys::BouncingBall, x)
    Wn = wn(sys, x)
    return Wn'*x[3:4]
end

function vtang(sys::BouncingBall, x)
    Wt = wt(sys, x)
    return Wt'*x[3:4]
end

function massMatrix(sys::BouncingBall, x)
    return [sys.m 0.0;
            0.0 sys.m]
end

function genForces(sys::BouncingBall, x)
    return [0.0, -sys.m*sys.g]
end

function wn(sys::BouncingBall, x)
    return permutedims(hcat([0.0, 1.0]...))
end

function wt(sys::BouncingBall, x)
    return permutedims(hcat([1.0, 0.0]...))
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
    plot(t, getindex.(Λn, 1))
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{n1} [N]$", fontsize=15)

    subplot(2, 1, 2)
    plot(t, getindex.(Λt, 1))
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{t1} [N]$", fontsize=15)
end
