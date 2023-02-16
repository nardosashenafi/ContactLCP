export RimlessWheel

const m1           = 2.0f0    
const m2           = 6.0f0
const l1           = 0.3f0        #wheel
const l2           = 0.06f0          #torso COM
const mt           = m1 + m2
const I1           = 0.0885f0
const I2           = m2*l2^2.0f0/3.0f0
const g            = 9.81f0
const k            = 10
const α            = Float32(360.0/k/2.0 * pi/180.0)
const γ            = Float32(0.0*pi/180.0)
const ls           = 5.0f0              #length of the slope (runway)
const ϵn_const     = 0.0f0*ones(Float32, k)
const ϵt_const     = 0.0f0*ones(Float32, k)
const μ_const      = 0.6f0*ones(Float32, k)
const gThreshold   = 0.001f0

struct RimlessWheel{}  
end

function initialState(θ0, θ0dot, ϕ0, ϕ0dot; γi = 0.0f0)
    @assert pi-α <= θ0 <= pi+α "Give an initial spoke angle for the spoke in contact. This helps calculate (x, y) correctly and sets the rimless wheel in contact with the surface. Pick initial θ0 such that pi-α <= θ0 <= pi+α"
    
    xs = 0.0f0
    ys = 0.0f0
    x = xs + l1*sin(θ0)
    y = ys - l1*cos(θ0) + 0.0001    #give it a little gap so it does not start with penetrating the groud
    ϕ = ϕ0 
    θ = θ0 
    xdot = l1*cos(θ0) * θ0dot
    ydot = l1*sin(θ0) * θ0dot
    ϕdot = ϕ0dot 
    θdot = θ0dot

    return [x, y, ϕ, θ, xdot, ydot, ϕdot, θdot]
end

#returns the attributes needed to model contact
function (sys::RimlessWheel)(z, param::Vector{T}; ϵn=ϵn_const, ϵt=ϵt_const, μ=μ_const, gThreshold=gThreshold, expert=false) where {T<:Real}
    q, v = parseStates(z)
    x, y, ϕ, θ = q
    gn, Wn, Wt = gap(z)
    γn  = vnormal(Wn, v)
    γt  = vtang(Wt, v)
    M   = massMatrix(ϕ)
    h   = genForces(z, param; expert=expert)

    return gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold
end

function (sys::RimlessWheel)(z, sysParam, controlParam::Vector{T}; ϵn=ϵn_const, ϵt=ϵt_const, μ=μ_const, gThreshold=gThreshold, expert=false) where {T<:Real}
    q, v = parseStates(z)
    x, y, ϕ, θ = q
    gn, Wn, Wt = gap(z)
    γn  = vnormal(Wn, v)
    γt  = vtang(Wt, v)
    M   = massMatrix(ϕ)
    h   = genForces(z, controlParam; expert=expert)

    return gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold
end

#parses out the poses and velocities from the state vector
function (sys::RimlessWheel)(x::Vector{T}) where {T<:Real}
    return parseStates(x)
end

function parseStates(x::Vector{T}) where {T<:Real}
    q = x[1:4]      #[x, y, ϕ, θ], ϕ is torso angle and θ is angle of the spoke in contact. If there are two in contact, it can be either one. Gap function checks either anyway
    v = x[5:8]      #[xdot, ydot, ϕdot, θdot]
    return q, v
end

function gap(z, γi)     #allows you to change the slope but contact force is not two force member here. 
    
    k_range         = range(0, stop=k-1, step=1)
    xhip, yhip, θ   = (z[1], z[2], z[4])
    gn              = zeros(k)

    for ki in k_range 
        xs, ys          = (xhip - l1*sin(θ + 2*α*ki), yhip + l1*cos(θ + 2*α*ki) )
        x_intersection  = (xs*yhip - ys*xhip)/(tan(γi)*(xs - xhip) - (ys - yhip))
        y_intersection  = tan(γi)*x_intersection            #intersection between the spokes and the incline
        gn[ki+1]        = ys - y_intersection
    end

    return gn
end

function gap(z::Vector{T}) where {T<:Real}     #assumes flat surface
    
    k_range         = range(0, stop=k-1, step=1)
    xhip, yhip, θ   = (z[1], z[2], z[4])
    gn              = zeros(T, k)
    Wn              = zeros(T, 4, k)
    Wt              = zeros(T, 4, k)
    for ki in k_range 

        θi          = θ + 2.0f0*α*ki
        xs, ys      = (xhip - l1*sin(θ + 2.0f0*α*ki), yhip + l1*cos(θi) )
        x_floor     = (ys*xhip - xs*yhip) / (ys-yhip)
        gn[ki+1]    = sign(ys)*sqrt((xs - x_floor)^2 + ys^2)
        Wn[:, ki+1] = 1.0f0/(gn[ki+1]+1e-5) .* [(xs-x_floor) ys 0.0f0 l1*cos(θi)*(x_floor-xs)-ys*l1*sin(θi)]
        Wt[:, ki+1] = [1.0f0 0.0f0 0.f0 -l1*cos(θi)]

    end
    return gn, Wn, Wt
end

function vnormal(Wn, v)
    return Wn'*v
end

function vtang(Wt, v)
    return Wt'*v
end

function massMatrix(ϕ)
    
    M = [mt 0.0f0 m2*l2*cos(ϕ) 0.0f0;
        0.0f0 mt m2*l2*sin(ϕ) 0.0f0;
        m2*l2*cos(ϕ) m2*l2*sin(ϕ) I2+m2*l2^2 0.0f0;
        0.0f0 0.0f0 0.0f0 I1]

    return M
end

inputLayer(x) = [cos(x[3]), sin(x[3]), cos(x[4]), sin(x[4]), x[7], x[8]]

function control(z, θp::Vector{T}; expert=false) where {T<:Real}
    q, v = parseStates(z)

    if expert #working expert controller
        @assert length(θp) == 2
        return -θp[1]*(q[3]-0.38f0) - θp[2]*v[3]
    else
        # u = -θp[1]*(q[3]-0.38f0) - θp[2]*v[3]
        # @assert length(θp) == DiffEqFlux.paramlength(unn) 
        # u = unn(inputLayer(z), θp)[1]
        @assert length(θp) == DiffEqFlux.paramlength(Hd) + N
        u = MLBasedESC.controller(npbc, inputLayer(z), θp)
        
        return clamp(u, -satu, satu)
    end
end

function genForces(z, param::Vector{T}; expert=false) where {T<:Real}

    q, v = parseStates(z)
    x, y, ϕ, θ = q
    xdot, ydot, ϕdot, θdot = v
    #h = Bu - C qdot - G
    B = [0.0f0, 0.0f0, 1.0f0, -1.0f0]

    C = [0.0f0 0.0f0 -m2*l2*sin(ϕ)*ϕdot 0.0f0;
        0.0f0 0.0f0 m2*l2*cos(ϕ)*ϕdot 0.0f0;
        -m2*l2*sin(ϕ)*ϕdot m2*l2*cos(ϕ)*ϕdot 0.0f0 0.0f0;
        0.0f0 0.0f0 0.0f0 0.0f0]

    G = [mt*g*sin(γ), 
        mt*g*cos(γ), 
        m2*g*l2*sin(ϕ + γ),
        0.0f0]

    return B*control(z, param; expert=expert) - C*v - G
end

function genForces(z, γi, param::Vector{T}; expert=false) where {T<:Real}

    q, v = parseStates(z)
    x, y, ϕ, θ = q
    xdot, ydot, ϕdot, θdot = v
    #h = Bu - C qdot - G
    B = [0.0f0, 0.0f0, 1.0f0, -1.0f0]

    C = [0.0f0 0.0f0 -m2*l2*sin(ϕ)*ϕdot 0.0f0;
        0.0f0 0.0f0 m2*l2*cos(ϕ)*ϕdot 0.0f0;
        -m2*l2*sin(ϕ)*ϕdot m2*l2*cos(ϕ)*ϕdot 0.0f0 0.0f0;
        0.0f0 0.0f0 0.0f0 0.0f0]

    G = [-mt*g*sin(γi), 
        mt*g*cos(γi), 
        m2*g*l2*sin(ϕ - γi),
        0.0f0]

    return B*control(z, param; expert=expert) - C*v - G
end

function impactMap(ϕ)

    det = I1*I2 + I1*m2*l2^2 + I2*mt*l1^2 + 
            m2*l1^2*l2^2*(m1 + m2*(sin(α - ϕ))^2)

    ξ1 = 1/det * ((I1*I2 + I1*m2*l2^2) + 
            (I2*mt*l1^2 + m2*l1^2*l2^2*(m1 +0.5* m2))*cos(2α) -
            1/2*m2^2*l1^2*l2^2*cos(2ϕ))

    ξ2 = 1/det *(m2*l1*l2*(I1*(cos(α - ϕ) - cos(α + ϕ)) + 
            mt*l1^2*(cos(2α)*cos(α - ϕ) - cos(α + ϕ))) )

    return [ξ1 0;
            ξ2 1]

end

function comparePostImpact(lcp, x, param; Δt = 0.001f0)
    
    ϕ = x[3]
    θdot, ϕdot = [x[8], x[7]]
    postImpactMap = impactMap(ϕ) * [θdot, ϕdot]

    gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold = sysAttributes(lcp, x, param)
    λn, λt, λR  = solveLcp(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold, x; Δt=Δt)

    qM, uA = sys(x)
    uE = M\((Wn - Wt*diagm(0 => lcp.μ))*λn + Wt*λR + h*Δt) + uA
    θdotlcp, ϕdotLcp  = [uE[4], uE[3]]

    println("[θdot-, ϕdot-] = ", [θdot, ϕdot])
    println("Impact map = ", postImpactMap)
    println("LCP = ", [θdotlcp, ϕdotLcp])

end

function createAnimateObject(x, y, ϕ, θ; spokeGap=0.2, k=k, α=α, γ=γ)

    vspokes1 = vis[:spokes1]

    for ki in range(0, stop=k-1, step=1)
        vki = vspokes1[Symbol("spoke" * String("$ki"))]

        setobject!(vki, MeshObject(
            Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, l1), 0.015),
            MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.7))))
        settransform!(vki, Translation(0.0, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
    end

    vspokes2 = vis[:spokes2]

    for ki in range(0, stop=k-1, step=1)
        vki = vspokes2[Symbol("spoke" * String("$ki"))]

        setobject!(vki, MeshObject(
            Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, l1), 0.015),
            MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.25))))
        settransform!(vki, Translation(spokeGap, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
    end

    vtorso = vis[:torso]
    setobject!(vtorso[:link], MeshObject(
        HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(spokeGap, 0.02, l2+0.1)),
        MeshLambertMaterial(color=RGBA{Float32}(1.0, 0.0, 0.0, 1.0))))
    settransform!(vtorso[:link], Translation(0.0, x+0.01, y) ∘ LinearMap(RotX(ϕ+pi)))
   
    vrunway = vis[:runway]
    setobject!(vrunway[:runway], MeshObject(
        Rect(Vec(0.0, 0.0, 0.0), Vec(spokeGap, ls, 0.01)),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.7))))
    settransform!(vrunway[:runway], Translation(0.0, x, -ls*sin(γ)) ∘ LinearMap(RotX(γ)))

    return vspokes1, vspokes2, vtorso
end

function animate(Z; γ=γ)

    x0, y0, ϕ0, θ0 = Z[1][1:4]
    spokeGap = 0.1
    vspokes1, vspokes2, vtorso = createAnimateObject(x0, y0, ϕ0, θ0; γ=γ, spokeGap=spokeGap)
    for z in Z[1:30:end]
        x, y, ϕ, θ = z[1:4]
        for ki in range(0, stop=k-1, step=1)
            vki = vspokes1[Symbol("spoke" * String("$ki"))]
            settransform!(vki, Translation(0.0, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
        end
        for ki in range(0, stop=k-1, step=1)
            vki = vspokes2[Symbol("spoke" * String("$ki"))]
            settransform!(vki, Translation(spokeGap, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
        end
        settransform!(vtorso[:link], Translation(0.0, x+0.01, y) ∘ LinearMap(RotX(ϕ+pi)))
        sleep(0.005)
    end

end

function startAnimator()
    window = Window()
    vis = Visualizer()
    open(vis, window)

    return vis
end

function animate_random(lcp, param; totalTimeStep = 5000)
    x0 = Float32.(initialState(rand(pi-α:0.01:pi+α), 
                                        #  rand(-3.0:0.05:-0.5),
                                        0.0, 
                                         0.0, 
                                         0.0) )

    X, _, _, _ = fulltimestep(lcp, x0, param; Δt = 0.001f0, totalTimeStep = totalTimeStep);                            
    animate(X)
end

function plots(Z)
    fig1 = plt.figure(1)
    plots(Z, fig1)

end

function plots(Z, fig1)
    PyPlot.figure(1)
    fig1.clf()
    subplot(2, 2, 1)
    plot(getindex.(Z, 3), getindex.(Z, 7))
    ylabel(L"\dot{\phi} [rad/s]", fontsize=15)
    subplot(2, 2, 2)
    θ, impactIndex = spokeInContact(Z)
    plot(θ, getindex.(Z, 8))
    scatter(θ[end], Z[end][8])
    ylabel(L"\dot{\theta} [rad/s]", fontsize=15)
    subplot(2, 2, 3)
    plot(getindex.(Z, 5))
    ylabel("vx [m/s]", fontsize=15)
    println("Average hip speed = ", mean(getindex.(Z, 5)))
end

function plots(Z, t, Λn, Λt)
    fig1 = plt.figure(1)
    fig2 = plt.figure(2)
    fig3 = plt.figure(3)
    plots(Z, t, Λn, Λt, fig1, fig2, fig3)
end

function plots(Z, t, Λn, Λt, fig1, fig2, fig3)
    PyPlot.figure(1)
    fig1.clf()
    subplot(2, 3, 1)
    plot(t, getindex.(Z, 1))
    ylabel("x [m]", fontsize=15)
    subplot(2, 3, 2)
    plot(t, getindex.(Z, 2))
    ylabel("y [m]", fontsize=15)
    subplot(2, 3, 3)
    plot(t, getindex.(Z, 5))
    ylabel("vx [m/s]", fontsize=15)
    subplot(2, 3, 4)
    plot(t, getindex.(Z, 6))
    ylabel("vy [m/s]", fontsize=15)
    subplot(2, 3, 5)
    plot(getindex.(Z, 3), getindex.(Z, 7))
    ylabel(L"\dot{\phi} [rad/s]", fontsize=15)
    subplot(2, 3, 6)
    plot(getindex.(Z, 4), getindex.(Z, 8))
    ylabel(L"\dot{\theta} [rad/s]", fontsize=15)

    PyPlot.figure(2)
    fig2.clf()
    subplot(2, 1, 1)
    for i in 1:length(Λn[1])
        plot(t, getindex.(Λn, i))
    end
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{n1} [N]$", fontsize=15)

    subplot(2, 1, 2)
    for i in 1:length(Λt[1])
        plot(t, getindex.(Λt, i))
    end
    ticklabel_format(axis="y", style="sci",scilimits=(0,0))
    ylabel(L"$\lambda_{t1} [N]$", fontsize=15)

    PyPlot.figure(3)
    fig3.clf()
    θ, impactIndex = spokeInContact(Z)
    plot(θ, getindex.(Z, 8))
end

function spokeInContact(Z; gThreshold=gThreshold, k=k, α=α, Δt=0.0005f0)

    θ = getindex.(Z, 4)
    θ_new = deepcopy(θ)
    θdot = getindex.(Z, 8)
    i     = 1
    q, v = parseStates(Z[1])
    znew = vcat(q + v*Δt/2.0f0, v)
    gn = gap(znew)
    Wn = wn(q[4])
    γn = vnormal(Wn, v)
    contactIndex, _ = checkContact(znew, gn, gThreshold, k)
    ki_prev = contactIndex[1]
    impactInd = Vector{Int32}(undef, length(Z))

    for z in Z
        q, v = parseStates(z)
        znew = vcat(q + v*Δt/2.0f0, v)
        gn = gap(znew)
        Wn = wn(q[4])
        γn = vnormal(Wn, v)
        contactIndex, _ = checkContact(znew, gn, gThreshold, k)
        
        if !isempty(contactIndex) 
            if !isnothing(findfirst(γni -> γni < 0.0, γn[contactIndex]))
                ki = contactIndex[findfirst(γni -> γni < 0.0, γn[contactIndex])]   #the small and decreasing gap
                ki != ki_prev ? impactInd[i] = 1 : impactInd[i] = 0
                if 0 <= abs(ki - ki_prev) < k-1
                    δk = (ki - ki_prev) 
                elseif abs(ki - ki_prev) == k-1
                    δk = sign(ki - ki_prev)
                end
                θ_new[i+1:end] .+= 2*α*δk
                ki_prev = ki
            end
        end
        i += 1
    end

    return θ_new, impactInd
end

