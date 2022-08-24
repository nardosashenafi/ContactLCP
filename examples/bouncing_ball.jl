struct BouncingBall{T}
    m               ::T 
    r               ::T
    g               ::T 
    ϵn              ::Vector{T}
    ϵt              ::Vector{T}
    μ               ::Vector{T}
    x0              ::Vector{T}
    contactIndex   ::Vector{T}

    function BouncingBall(T)

        m               = T(0.2)
        r               = T(0.1)
        g               = T(9.81)
        ϵn              = T.(0.9*ones(1))
        ϵt              = T.(-0.5*ones(1))
        μ               = T.(0.2*ones(1))
        x0              = T.([0.0, 0.5, 0.1, -0.1])     #xpos, ypos, xdot, ydot
        contactIndex   = zeros(T, 1)

        new{T}(m, r, g, ϵn, ϵt, μ, x0, contactIndex)
    end

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

function massMatrix(sys, x)
    return [sys.m 0.0;
            0.0 sys.m]
end

function genForces(sys, x)
    return [0.0, -sys.m*sys.g]
end

function wn(sys, x)
    return permutedims(hcat([0.0, 1.0]...))
end

function wt(sys, x)
    return permutedims(hcat([1.0, 0.0]...))
end