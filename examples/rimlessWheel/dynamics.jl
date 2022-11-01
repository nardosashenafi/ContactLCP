export RimlessWheel

const m1           = 2.32f0    
const m2           = 4.194f0 
const I1           = 0.0784160f0
const I2           = 0.0380256f0
const mt           = m1 + m2
const l1           = 0.26f0        #wheel
const l2           = 0.05f0          #torso
const g            = 9.81f0
const k            = 10
const α            = Float32(360.0/k/2.0 * pi/180.0)
const γ            = Float32(0.0*pi/180.0)
const ϵn_const     = 0.0f0*ones(Float32, k)
const ϵt_const     = 0.0f0*ones(Float32, k)
const μ_const      = 0.6f0*ones(Float32, k)
const gThreshold   = 0.001f0

struct RimlessWheel{}  
end

function initialState(θ0, θ0dot, ϕ0, ϕ0dot)
    @assert pi-α <= θ0 <= pi+α "Give an initial spoke angle for the spoke in contact. This will help set the rimless wheel in contact with the surface"
    
    x = l1*sin(pi-θ0)
    y = l1*cos(pi-θ0)
    ϕ = ϕ0 
    θ = θ0 
    xdot = -l1*cos(pi-θ0) * θ0dot
    ydot = l1*sin(pi-θ0) * θ0dot
    ϕdot = ϕ0dot 
    θdot = θ0dot

    return [x, y, ϕ, θ, xdot, ydot, ϕdot, θdot]
end

#returns the attributes need to model contact
function (sys::RimlessWheel)(x, θ::Vector{T}; ϵn=ϵn_const, ϵt=ϵt_const, μ=μ_const, gThreshold=gThreshold, expert=false) where {T<:Real}
    gn  = gap(x)  
    γn  = vnormal(x)
    γt  = vtang(x)
    M   = massMatrix(x)
    h   = genForces(x, θ; expert=expert)
    Wn  = wn(x)
    Wt  = wt(x)

    return gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold
end

#parses out the poses and velocities from the state vector
function (sys::RimlessWheel)(x::Vector{T}) where {T<:Real}
    return parseStates(x)
end

function parseStates(x::Vector{T}) where {T<:Real}
    q = x[1:4]      #[x, y, ϕ, θ]
    u = x[5:8]      #[xdot, ydot, ϕdot, θdot]
    return q, u
end

function gap(z)
    
    k_range = range(0, stop=k-1, step=1)
    y, θ = (z[2], z[4])
    return y .+ l1*cos.(θ .+ 2*α*k_range) 
end

function vnormal(z)
    Wn = wn(z)
    _, v = parseStates(z)
    return Wn'*v
end

function vtang(z)
    Wt = wt(z)
    _, v = parseStates(z)
    return Wt'*v
end

function massMatrix(z)
    q, v = parseStates(z)
    x, y, ϕ, θ = q
    
    M = [mt 0.0 m2*l2*cos(ϕ) 0.0;
        0.0 mt m2*l2*sin(ϕ) 0.0;
        m2*l2*cos(ϕ) m2*l2*sin(ϕ) I2+m2*l2^2 0.0;
        0.0 0.0 0.0 I1]

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

    G = [-mt*g*sin(γ), 
        mt*g*cos(γ), 
        m2*g*l2*sin(ϕ - γ),
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

function comparePostImpact(lcp, x, param; Δt = 0.001)
    
    ϕ = x[3]
    θdot, ϕdot = [x[8], x[7]]
    postImpactMap = impactMap(ϕ) * [θdot, ϕdot]

    gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold = sysAttributes(lcp, x, param)
    λn, λt, λR  = solveLcp(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold; Δt=Δt)

    qM, uA = sys(x)
    uE = M\((Wn - Wt*diagm(0 => lcp.μ))*λn + Wt*λR + h*Δt) + uA
    θdotlcp, ϕdotLcp  = [uE[4], uE[3]]

    println("[θdot-, ϕdot-] = ", [θdot, ϕdot])
    println("Impact map = ", postImpactMap)
    println("LCP = ", [θdotlcp, ϕdotLcp])

end

function wn(z::Vector{T}) where {T<:Real}

    θ        = z[4]
    k_range  = range(0, stop=k-1, step=1)
    Wn       = zeros(T, 4, k)
    Wn[2, :] = ones(T, k)
    Wn[4, :] = -l1 .* sin.(θ .+ 2.0*α.*k_range)

   return Wn 

end

function wt(z::Vector{T}) where {T<:Real}

    θ        = z[4]
    k_range  = range(0, stop=k-1, step=1)
    Wt       = zeros(T, 4, k)
    Wt[1, :] = ones(T, k)
    Wt[4, :] = -l1 .* cos.(θ .+ 2.0*α.*k_range)

    return Wt
end

function createAnimateObject(x, y, ϕ, θ; k=k, α=α)
    vspokes = vis[:spokes]

    for ki in range(0, stop=k-1, step=1)
        vki = vspokes[Symbol("spoke" * String("$ki"))]

        setobject!(vki, MeshObject(
            Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, l1), 0.015),
            MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.25))))
        settransform!(vki, Translation(0.0, x, y) ∘ LinearMap(RotX(θ+2*α*ki)))
    end

    vtorso = vis[:torso]
    setobject!(vtorso[:link], MeshObject(
        Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, -l2), 0.005),
        MeshLambertMaterial(color=RGBA{Float32}(1.0, 0.0, 0.0, 1.0))))
    settransform!(vtorso[:link], Translation(0.0, x, y) ∘ LinearMap(RotX(ϕ)))

    setobject!(vtorso[:bob], MeshObject(
        HyperSphere(Point(0.0, 0.0, -l2), 0.015),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 1.0, 1.0))))
    settransform!(vtorso[:bob], Translation(0.0, x, y) ∘ LinearMap(RotX(ϕ)))
   
    return vspokes, vtorso
end

function animate(Z)

    x0, y0, ϕ0, θ0 = Z[1][1:4]
    vspokes, vtorso = createAnimateObject(x0, y0, ϕ0, θ0)
    for z in Z[1:50:end]
        x, y, ϕ, θ = z[1:4]
        for ki in range(0, stop=k-1, step=1)
            vki = vspokes[Symbol("spoke" * String("$ki"))]
            settransform!(vki, Translation(0.0, x, y) ∘ LinearMap(RotX(θ+2*α*ki)))
        end
        settransform!(vtorso[:link], Translation(0.0, x, y) ∘ LinearMap(RotX(ϕ)))
        settransform!(vtorso[:bob], Translation(0.0, x, y) ∘ LinearMap(RotX(ϕ)))
        sleep(0.04)
    end

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
    θ = spokeInContact(Z)
    plot(θ, getindex.(Z, 8))
    ylabel(L"\dot{\theta} [rad/s]", fontsize=15)
    subplot(2, 2, 3)
    plot(getindex.(Z, 5))
    ylabel("vx [m/s]", fontsize=15)
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
    θ = spokeInContact(Z)
    plot(θ, getindex.(Z, 8))
end

# function spokeInContact(Z; gThreshold=gThreshold, k=k, α=α)

#     θ = getindex.(Z, 4)
#     θ_new = deepcopy(θ)
#     i     = 1

#     for z in Z
#         gn = gap(z)
#         contactIndex, _ = checkContact(gn, gThreshold, k)

#         if any(contactIndex .> 0.0)
#             ki = findfirst(x -> x == 1.0, contactIndex) - 1
#             θ_new[i:end] = θ[i:end] .+ 2*α*ki
#         end
#         i += 1
#     end

#     return θ_new
# end


function spokeInContact(Z)
    θ_new = getindex.(Z, 4)
    θdot = getindex.(Z, 8)

    for i in 2:length(θ_new)
        #if velocity jumps, wrap. Checking θ causes the angle to jump when the velocity has not
        if (abs(θdot[i] - θdot[i-1]) > 0.1) && 
            θdot[i] < 0.0  &&
            abs(round(abs.(θ_new[i]/α)) - abs.(θ_new[i]/α)) < 0.03       #instead check if new contact occurs
            
            θ_new[i:end] .= θ_new[i:end] .+ 2*α

        elseif (abs(θdot[i] - θdot[i-1]) > 0.1) &&
            θdot[i] > 0.0 && 
            abs(round(abs.(θ_new[i]/α)) - abs.(θ_new[i]/α)) < 0.03      
            
            θ_new[i:end] .= θ_new[i:end] .- 2*α
        end
    end
    return θ_new
end
