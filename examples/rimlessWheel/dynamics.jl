export RimlessWheel

const m1           = 1.13f0    
const m2           = 3.385f0
const l1           = 0.3f0        #wheel
const l2           = 0.06f0          #torso COM
const mt           = m1 + m2
const I1           = 0.0885f0/2.0f0
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

function initialStateWithBumps(θ0::T, θ0dot, ϕ0, ϕ0dot, rmax::T; k=k, α=α) where {T<:Real}
    @assert pi-α <= θ0 <= pi+α "Give an initial spoke angle for the spoke in contact. This helps calculate (x, y) correctly and sets the rimless wheel in contact with the surface. Pick initial θ0 such that pi-α <= θ0 <= pi+α"
    
    r = rmax*rand(T, k)
    x = l1*sin(pi-θ0)
    y = r[1] + l1*cos(pi-θ0) 
    ϕ = ϕ0 
    θ = θ0 
    xdot = -l1*cos(pi-θ0) * θ0dot
    ydot = l1*sin(pi-θ0) * θ0dot
    ϕdot = ϕ0dot 
    θdot = θ0dot

    #make sure the next spoke is not penetrating into the bump
    r[2] = r[1]*0.5
    r[10] = r[1]*0.5
    return [x, y, ϕ, θ, xdot, ydot, ϕdot, θdot], r
end

function initialStateWithBumps(θ0::T, θ0dot, ϕ0, ϕ0dot, rvec::Vector{T}; k=k, α=α) where {T<:Real}
    @assert pi-α <= θ0 <= pi+α "Give an initial spoke angle for the spoke in contact. This helps calculate (x, y) correctly and sets the rimless wheel in contact with the surface. Pick initial θ0 such that pi-α <= θ0 <= pi+α"
    
    x = l1*sin(pi-θ0)
    y = rvec[1] + l1*cos(pi-θ0) 
    ϕ = ϕ0 
    θ = θ0 
    xdot = -l1*cos(pi-θ0) * θ0dot
    ydot = l1*sin(pi-θ0) * θ0dot
    ϕdot = ϕ0dot 
    θdot = θ0dot

    return [x, y, ϕ, θ, xdot, ydot, ϕdot, θdot], rvec
end

"""
    sys(state, controllerParameters)

Return system specific parameters such as gap function to model contact.

This call method is what relays the system information to the LCP 
"""

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

""" 
    sys(state, systemParameters, controllerParameters)

Return system specific parameters such as gap function to model contact.

Pass system parameters to introduce random bumps to each spoke. "systemParameters" contains [[x1, x2...], [y1, y2,...], [height1, height2...]] of each bump for each spoke

"""
function (sys::RimlessWheel)(z, sysParam, controlParam::Vector{T}; ϵn=ϵn_const, ϵt=ϵt_const, μ=μ_const, gThreshold=gThreshold, expert=false) where {T<:Real}
    q, v = parseStates(z)
    x, y, ϕ, θ = q
    gn, Wn, Wt = gap(z, sysParam)
    γn  = vnormal(Wn, v)
    γt  = vtang(Wt, v)
    M   = massMatrix(ϕ)
    h   = genForces(z, controlParam; expert=expert)

    return gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold
end

""" 
    sys(state)
Parse out the poses and velocities from the state vector. Used by LCP

"""

function (sys::RimlessWheel)(x::Vector{T}) where {T<:Real}
    return parseStates(x)
end

"""
    createAnimateObject(state, bumpHeight)

Return bump location [[x1, x2...], [y1, y2,...], [height1, height2...]] for each spoke

For each spoke, there is a designated bump of height r. The height remains fixed throughout a trajectory.
This function places the bumps at the appriopriate location on the ground so the spokes can strike the designated bumps.

"""
function createUnevenTerrain(z::Vector{T}, r) where {T<:Real}
    x, y, θ = (z[1], z[2], z[4])
    ki = spokeNearGround(gap(z)[1]) 
    xs, ys   = (x - l1*sin(θ + 2.0*α*(ki-1)), y + l1*cos(θ + 2.0*α*(ki-1)) )

    xb = zeros(T, k)
    yb = zeros(T, k)

    xb[ki] = xs 
    kj_prev = 0

    for kj in vcat(range(ki+1, k, step=1), range(1, ki-1, step=1))
        kj-1 < 1 ? kj_prev = 10 : kj_prev = kj-1
        l̄      = sqrt((l1)^2 + (l1)^2 - 2*(l1)*(l1)*cos(2α))
        xb[kj] = xb[kj_prev] + l̄
    end

    return [xb, yb, r]
end

function parseStates(x::Vector{T}) where {T<:Real}
    q = x[1:4]      #[x, y, ϕ, θ], ϕ is torso angle and θ is angle of the spoke in contact. If there are two in contact, it can be either one. Gap function checks either anyway
    v = x[5:8]      #[xdot, ydot, ϕdot, θdot]
    return q, v
end

function gap(z::Vector{T}) where {T<:Real}      
    # gap to a level floor; the gap is the vertical distance between end of spoke and ground
    k_range = range(0, stop=k-1, step=1)
    y, θ     = (z[2], z[4])
    gn       = y .+ l1*cos.(θ .+ 2*α*k_range) 
    Wn       = zeros(T, 4, k)
    Wn[2, :] = ones(T, k)
    Wn[4, :] = -l1 .* sin.(θ .+ 2.0*α.*k_range)

    Wt       = zeros(T, 4, k)
    Wt[1, :] = ones(T, k)
    Wt[4, :] = -l1 .* cos.(θ .+ 2.0*α.*k_range)

    return gn, Wn, Wt
end

function hangingSpoke(z, gn; gThreshold=gThreshold, k=k)
    hangingSpoke   = 0
    contactIndex,_ = checkContact(z, gn, gThreshold, k)

    if z[8] <= 0.0
        hangingSpoke = maximum(contactIndex) + 1
        hangingSpoke > k ? hangingSpoke=1 : nothing 
    else
        hangingSpoke = minimum(contactIndex) - 1
        hangingSpoke < 1 ? hangingSpoke=k : nothing 
    end
end

"""
    gap(state, bumpParameters)
    
Compute gap between edge of spokes and corresponding bump.  
See also [`createUnevenTerrain`](@ref) on how to construct bumps  

Constructs constant bump for each spoke; the gap is measured along the length of the spoke

"""

function gap(z::Vector{T}, sysParam)  where {T<:Real}
    xb, yb, r = sysParam
    gn       = zeros(T, k)
    Wn       = zeros(T, 4, k)
    Wt       = zeros(T, 4, k)

    for ki in range(0, stop=k-1, step=1)     
        gki, wnki, wtki = gapSpokeToObstacle(z, ki, xb[ki+1], yb[ki+1], r[ki+1])
        gn[ki+1] = gki
        Wn[:, ki+1] = wnki
        Wt[:, ki+1] = wtki
    end

    return gn, Wn, Wt
end

function gapSpokeToObstacle(z::Vector{T}, ki, xbi, ybi, ri) where {T<:Real}     
    #takes gap along the length of the spoke
    xhip, yhip, θ   = (z[1], z[2], z[4])
    δ               = 1e-5          #for divide by zeros

    θi       = θ + 2.0f0*α*ki
    xs, ys   = (xhip - l1*sin(θi), yhip + l1*cos(θi) )

    ####find point of contact on the bump 
    # θb      = atan(ys-ybi, xs-xbi)
    # x_floor = xbi + ri*cos(θb)
    # y_floor = ybi + ri*sin(θb)
    x_floor = xbi
    y_floor = ri
    #########compute gap and its derivatives
    if  xbi - 0.05 <= xs <= xbi + 0.05
        # gn = sign(ys-y_floor)*sqrt((xs - x_floor)^2 + (ys-y_floor)^2)
        # Wn = 1.0f0/(gn+δ) .* [(xs-x_floor) (ys-y_floor) 0.0f0 l1*cos(θi)*(x_floor-xs)-(ys-y_floor)*l1*sin(θi)]
        # Wt = [ys-y_floor x_floor-xs 0.0f0 -l1*cos(θi)*(ys-y_floor)-l1*sin(θi)*(x_floor-xs)]

        gn = ys - y_floor
        Wn = [0.0f0 1.0f0 0.0f0 -l1*sin(θ + 2.0*α*ki)]
        Wt = [1.0f0 0.0 0.0 -l1*cos(θ + 2.0*α*ki)]

        return gn, Wn, Wt
    else 
        gn = ys
        Wn = [0.0f0 1.0f0 0.0f0 -l1*sin(θ + 2.0*α*ki)]
        Wt = [1.0f0 0.0 0.0 -l1*cos(θ + 2.0*α*ki)]
        return gn, Wn, Wt
    end
end

function spokeNearGround(gn)
    _, spoke = findmin(abs.(gn))
    return spoke
end

function vnormal(Wn, v)
    return Wn'*v
end

function vtang(Wt, v)
    return Wt'*v
end

function massMatrix(ϕ)
    
    return [mt 0.0f0 m2*l2*cos(ϕ) 0.0f0;
            0.0f0 mt m2*l2*sin(ϕ) 0.0f0;
            m2*l2*cos(ϕ) m2*l2*sin(ϕ) I2+m2*l2^2 0.0f0;
            0.0f0 0.0f0 0.0f0 I1]

end

inputLayer(x) = [cos(x[3]), sin(x[3]), cos(x[4]), sin(x[4]), x[7], x[8]]

function control(z, θp::Vector{T}; expert=false) where {T<:Real}
    q, v = parseStates(z)

    if expert #working expert controller
        @assert length(θp) == 2
        return -θp[1]*(q[3]-0.65f0) - θp[2]*v[3]
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

function comparePostImpact(lcp, x; param=zeros(2), expert=true, Δt = 0.001f0)

    println("Make sure x is a point of impact")
    # x= [0.09163268650073744, 0.28589820307479014, -0.006866646029560244, 2.8311059731226207, 1.6253018883308197, -0.5162607562984872, 0.013881086255179611, -5.657496313791293]
    ϕ = x[3]
    θdot, ϕdot = [x[8], x[7]]
    postImpactMap = impactMap(ϕ) * [θdot, ϕdot]

    qA, uA  = lcp.sys(x)
    qM      = qA + 0.5f0*Δt*uA

    x_mid   = vcat(qM, uA)
    gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold = sysAttributes(lcp, x_mid, param;expert=expert)
    λn, λt, λR  = solveLcp(gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold, x_mid; Δt=Δt)

    uE = M\((Wn - Wt*diagm(0 => μ_const))*λn + Wt*λR + h*Δt) + uA
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

    # vspokes2 = vis[:spokes2]

    # for ki in range(0, stop=k-1, step=1)
    #     vki = vspokes2[Symbol("spoke" * String("$ki"))]

    #     setobject!(vki, MeshObject(
    #         Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, l1), 0.015),
    #         MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.25))))
    #     settransform!(vki, Translation(spokeGap, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
    # end

    vtorso = vis[:torso]
    setobject!(vtorso[:link], MeshObject(
        HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(spokeGap, 0.02, l2+0.1)),
        MeshLambertMaterial(color=RGBA{Float32}(1.0, 0.0, 0.0, 1.0))))
    settransform!(vtorso[:link], Translation(0.0, x+0.01, y) ∘ LinearMap(RotX(ϕ+pi)))
   
    # vrunway = vis[:runway]
    # setobject!(vrunway[:runway], MeshObject(
    #     Rect(Vec(0.0, 0.0, 0.0), Vec(spokeGap, ls, 0.01)),
    #     MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.7))))
    # settransform!(vrunway[:runway], Translation(0.0, x, -ls*sin(γ)) ∘ LinearMap(RotX(γ)))

    return vspokes1, vtorso #, vspokes2
end

function animate(Z; γ=γ)

    x0, y0, ϕ0, θ0 = Z[1][1:4]
    spokeGap = 0.1
    vspokes1, vtorso = createAnimateObject(x0, y0, ϕ0, θ0; γ=γ, spokeGap=spokeGap)
    # vspokes1, vtorso, vspokes2 = createAnimateObject(x0, y0, ϕ0, θ0; γ=γ, spokeGap=spokeGap)
    for z in Z[1:15:end]
        x, y, ϕ, θ = z[1:4]
        for ki in range(0, stop=k-1, step=1)
            vki = vspokes1[Symbol("spoke" * String("$ki"))]
            settransform!(vki, Translation(0.0, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
        end
        # for ki in range(0, stop=k-1, step=1)
        #     vki = vspokes2[Symbol("spoke" * String("$ki"))]
        #     settransform!(vki, Translation(spokeGap, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
        # end
        settransform!(vtorso[:link], Translation(0.0, x+0.01, y) ∘ LinearMap(RotX(ϕ+pi)))
        sleep(0.002)
    end
end

function createAnimateObject(x, y, ϕ, θ, sysParam; spokeGap=0.2, k=k, α=α, γ=γ)

    vspokes1 = vis[:spokes1]

    for ki in range(0, stop=k-1, step=1)
        vki = vspokes1[Symbol("spoke" * String("$ki"))]

        setobject!(vki, MeshObject(
            Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, l1), 0.015),
            MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.7))))
        settransform!(vki, Translation(0.0, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
    end

    # vspokes2 = vis[:spokes2]

    # for ki in range(0, stop=k-1, step=1)
    #     vki = vspokes2[Symbol("spoke" * String("$ki"))]

    #     setobject!(vki, MeshObject(
    #         Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, l1), 0.015),
    #         MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.25))))
    #     settransform!(vki, Translation(spokeGap, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
    # end

    vtorso = vis[:torso]
    setobject!(vtorso[:link], MeshObject(
        HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(spokeGap, 0.02, l2+0.1)),
        MeshLambertMaterial(color=RGBA{Float32}(1.0, 0.0, 0.0, 1.0))))
    settransform!(vtorso[:link], Translation(0.0, x+0.01, y) ∘ LinearMap(RotX(ϕ+pi)))
   
    xb, yb, r = sysParam
    vBumps = vis[:bumps]
    for ki in range(0, stop=k-1, step=1)
        vki = vBumps[Symbol("spoke" * String("$ki"))]
        setobject!(vki, MeshObject(
            HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(0.01, 0.1, r[ki+1])),
            MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 1.0, 0.75))))
        settransform!(vki, Translation(0.0, xb[ki+1]-0.05, yb[ki+1]) )
    end

    return vspokes1, vtorso, vBumps #, vspokes2
end

function animate(Z, sysParam; γ=γ)

    x0, y0, ϕ0, θ0 = Z[1][1:4]
    spokeGap = 0.1
    vspokes1, vtorso, vBumps = createAnimateObject(x0, y0, ϕ0, θ0, sysParam[1]; γ=γ, spokeGap=spokeGap)
    # vspokes1, vtorso, vspokes2 = createAnimateObject(x0, y0, ϕ0, θ0; γ=γ, spokeGap=spokeGap)
    for j in range(1, step=15, stop=length(Z))
        z = Z[j]
        x, y, ϕ, θ = z[1:4]
        for ki in range(0, stop=k-1, step=1)
            vki = vspokes1[Symbol("spoke" * String("$ki"))]
            settransform!(vki, Translation(0.0, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
        end
        # for ki in range(0, stop=k-1, step=1)
        #     vki = vspokes2[Symbol("spoke" * String("$ki"))]
        #     settransform!(vki, Translation(spokeGap, x, y) ∘ LinearMap(RotX(θ+2*α*ki+γ)))
        # end
        settransform!(vtorso[:link], Translation(0.0, x+0.01, y) ∘ LinearMap(RotX(ϕ+pi)))
        
        for ki in range(0, stop=k-1, step=1)
            vki = vBumps[Symbol("spoke" * String("$ki"))]
            settransform!(vki, Translation(0.0, sysParam[j][1][ki+1]-0.05, 0.0))
        end
        sleep(0.002)
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
    scatter(Z[end][3], Z[end][7])
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
    gn, Wn, _ = gap(znew)
    γn = vnormal(Wn, v)
    contactIndex, _ = checkContact(znew, gn, gThreshold, k)
    isempty(contactIndex) ? ki_prev = 1 : ki_prev=contactIndex[1]
    impactInd = Vector{Int32}(undef, length(Z))

    for z in Z
        q, v = parseStates(z)
        znew = vcat(q + v*Δt/2.0f0, v)
        gn, Wn, _ = gap(znew)
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

# function spokeInContact(Z; gThreshold=gThreshold, k=k, α=α, Δt=0.0005f0)

#     θ = getindex.(Z, 4)
#     θ_new = deepcopy(θ)
#     θdot = getindex.(Z, 8)
#     i     = 1
#     q, v = parseStates(Z[1])
#     znew = vcat(q + v*Δt/2.0f0, v)
#     gn, Wn, _ = gap(znew)
#     γn = vnormal(Wn, v)
#     contactIndex, _ = checkContact(znew, gn, gThreshold, k)
#     isempty(contactIndex) ? ki_prev = 1 : ki_prev=contactIndex[1]
#     impactInd = Vector{Int32}(undef, length(Z))

#     for z in Z
#         q, v = parseStates(z)
#         znew = vcat(q + v*Δt/2.0f0, v)
#         gn, Wn, _ = gap(znew)
#         γn = vnormal(Wn, v)
#         contactIndex, _ = checkContact(znew, gn, gThreshold, k)
        
#         if !isempty(contactIndex) 
#             if !isnothing(findfirst(γni -> γni < 0.0, γn[contactIndex]))
#                 ki = contactIndex[findfirst(γni -> γni < 0.0, γn[contactIndex])]   #the small and decreasing gap
#                 ki != ki_prev ? impactInd[i] = 1 : impactInd[i] = 0
#                 if 0 <= abs(ki - ki_prev) < k-1
#                     δk = (ki - ki_prev) 
#                 elseif abs(ki - ki_prev) == k-1
#                     δk = sign(ki - ki_prev)
#                 end
#                 θ_new[i+1:end] .+= 2*α*δk
#                 ki_prev = ki
#             end
#         end
#         i += 1
#     end

#     return θ_new, impactInd
# end

