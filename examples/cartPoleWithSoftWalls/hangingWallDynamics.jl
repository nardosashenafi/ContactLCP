using PyPlot 
using ControlSystems 
using MeshCat
using GeometryBasics
using CoordinateTransformations
using ColorTypes
using Blink
using Rotations

export CartPoleWithSoftWalls

const mc            = 1.0f0    
const mp            = 0.5f0 
const l             = 0.5f0        #wheel
const I1            = mp*l^2
const g             = 9.81f0
const d             = -0.75f0 
const wallThickness = 0.10f0
const D             = 1.5f0
const wallBottomEnd = 0.3f0
const wallTopEnd    = 0.8f0
const contactNum    = 4
const ϵn_const      = 0.5f0*ones(Float32, contactNum)
const ϵt_const      = 0.0f0*ones(Float32, contactNum)
const μ_const       = 0.0f0*ones(Float32, contactNum)
const gThreshold    = 0.001f0
const satu          = 4.0f0
const w             = 0.20f0

const leftWall_bottom   = [d, wallBottomEnd]
const leftWall_top      = [d, wallTopEnd]
const rightWall_bottom  = [D + d, wallBottomEnd]
const rightWall_top     = [D + d, wallTopEnd]

println("USING HANGING WALLS")
struct CartPoleWithSoftWalls{}  
end

function initialState(x1)
    # @assert pi-α <= θ0 <= pi+α "Give an initial spoke angle for the spoke in contact. This will help set the rimless wheel in contact with the surface"
    
    x1 = x1
    x2 = -10.0f0*pi/180f0

    x1dot = 0.0f0
    x2dot = 0.0f0

    return [x1, x2, x1dot, x2dot]
end

#returns the attributes needed to model contact
function (sys::CartPoleWithSoftWalls)(x, u::Vector{T}; ϵn=ϵn_const, ϵt=ϵt_const, μ=μ_const, gThreshold=gThreshold, expert=false) where {T<:Real}
    gn  = gap(x)  
    γn  = vnormal(x)
    γt  = vtang(x)
    M   = massMatrix(x)
    h   = genForces(x, u; expert=expert)
    Wn  = wn(x)
    Wt  = wt(x)

    return gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold
end

#parses out the poses and velocities from the state vector
function (sys::CartPoleWithSoftWalls)(x::Vector{T}) where {T<:Real}
    return parseStates(x)
end

function parseStates(x::Vector{T}) where {T<:Real}
    q = x[1:2]      #[x1, x2]
    u = x[3:4]      #[x1dot, x2dot]
    return q, u
end

function pendulumPos(z)
    q, _ = parseStates(z) 
    x, θ = q
    
    return [x-l*sin(θ), l*cos(θ)]
end

function gap(z)
    pendulum_xy = pendulumPos(z)

    if pendulum_xy[2] > leftWall_bottom[2] #pendulum in contact with left wall
        g1 = pendulum_xy[1] - d   
    else
        g1 = norm(pendulum_xy - leftWall_bottom, 2)
    end

    if pendulum_xy[2] > rightWall_bottom[2] #pendulum in contact with right wall
        g2 = (D + d) - pendulum_xy[1]
    else
        g2 = norm(pendulum_xy - rightWall_bottom, 2)
    end
    
    return [g1, g2, -g1 - wallThickness, -g2 - wallThickness] 
end

function wn(z::Vector{T}) where {T<:Real}

    q, _ = parseStates(z)
    x, θ = q
    pendulum_xy = pendulumPos(z)
    g = gap(z)

    if pendulum_xy[2] > leftWall_bottom[2] #pendulum in contact with left wall
        w1 = [1.0f0; -l*cos(q[2])]  
    else
        x̄l, ȳl = pendulum_xy - leftWall_bottom
        w1 = 2.0f0/g[1]*[x̄l; -x̄l*l*cos(θ) + ȳl*l*sin(θ)] 
    end

    if pendulum_xy[2] > rightWall_bottom[2] #pendulum in contact with right wall
        w2 = [-1.0f0; l*cos(q[2])]  
    else
        x̄r, ȳr = pendulum_xy - rightWall_bottom
        w2 = 2.0f0/g[2]*[x̄r; -x̄r*l*cos(θ) + ȳr*l*sin(θ)] 
    end
    
   return hcat(w1, w2, -w1, -w2)    #check contact against the wall on each face of the wall

end

function wt(z::Vector{T}) where {T<:Real}
    q, _ = parseStates(z)
    w1 = [0.0f0; -l*sin(q[2])]
    w2 = [0.0f0; -l*sin(q[2])]

    return hcat(w1, w2, w1, w2)
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
    q, u = parseStates(z)
    x1, x2 = q 
    
    M = [mc+mp -mp*l*cos(x2);
        -mp*l*cos(x2) mp*l^2+I1]

    return M
end

function lqrGains()
    M    = massMatrix([0.0f0, 0.0f0, 0.0f0, 0.0f0])
    Minv = inv(M)

    Ĝ = Minv* [0.0f0;
              mp*g*l]
    B̂ = Minv*[1.0f0, 0.0f0]

    A = [0.0f0 0.0f0 1.0f0 0.0f0;
        0.0f0  0.0f0 0.0f0 1.0f0;
        0.0f0  Ĝ[1]  0.0f0 0.0f0;
        0.0f0  Ĝ[2]  0.0f0 0.0f0]

    B = [0.0f0, 0.0f0, B̂[1], B̂[2]]

    C = Matrix{Float32}(LinearAlgebra.I, (4, 4)) 
    cartLinearized = ControlSystems.ss(A, B, C, 0.0f0)
    Q = 20.0f0*Matrix{Float32}(LinearAlgebra.I, (4, 4)) 
    R = 3.0f0
    ControlSystems.lqr(cartLinearized, Q, R)
end

function lqr(z)
    k = vec([ -2.58182  50.3793  -5.04367  12.0988])
    # k = zeros(4)
    return -k'*z
end

inputLayer(z) = [z[1], cos(z[2]), sin(z[2]), z[3], z[4]]
# inputLayer(z) = [cos(z[2]), sin(z[2]), z[3], z[4]]

function control(z, u::Vector{T}; expert=false, lqr_max = 10.0f0) where {T<:Real}
    q, v = parseStates(z)

    if expert #working expert controller
        return clamp(lqr(z), -lqr_max, lqr_max)
    else
        x1, x2 = q 
        x1dot, x2dot = v
        if ((1.0f0-cos(x2) <= (1.0f0-cosd(17.0))) && (abs(x2dot) <= 0.5f0))
            return clamp(lqr(z), -lqr_max, lqr_max)
        else
            return clamp(u[1], -satu, satu)
        end
    end
end

function genForces(z, u::Vector{T}; expert=false) where {T<:Real}

    q, v = parseStates(z)
    x1, x2 = q
    #h = Bu - C qdot - G
    B = [1.0f0, 0.0f0]

    C = [0.0f0 mp*l*sin(x2); 
        -mp*sin(x2)/2.0f0 0.0f0]

    G = [0.0f0;
        -mp*g*l*sin(x2)]

    return B*control(z, u; expert=expert) - C*v - G
end

function startAnimator()
    window = Window()
    vis = Visualizer()
    open(vis, window)

    return vis
end

function createAnimateObject(x1, x2)
    vcart = vis[:cart]

    setobject!(vcart, MeshObject(
        Rect(Vec(0.0, 0.0, 0.0), Vec(0.2, 0.2, 0.1)),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.25))))
    settransform!(vcart, Translation(-0.2, x1-w/2.0f0, 0.0))

    vpendulum = vis[:pendulum]
    setobject!(vpendulum[:link], MeshObject(
        Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, l), 0.005),
        MeshLambertMaterial(color=RGBA{Float32}(1.0, 0.0, 0.0, 1.0))))
    settransform!(vpendulum[:link], Translation(0.0, x1, 0.0) ∘ LinearMap(RotX(x2)))

    setobject!(vpendulum[:bob], MeshObject(
        HyperSphere(Point(0.0, 0.0, l), 0.015),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 1.0, 1.0))))
    settransform!(vpendulum[:bob], Translation(0.0, x1, 0.0) ∘ LinearMap(RotX(x2)))
   
    vwalls = vis[:walls]

    setobject!(vwalls[:left], MeshObject(
        Rect(Vec(0.0, 0.0, wallBottomEnd), Vec(0.0, wallThickness, wallTopEnd-wallBottomEnd)),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 1.0))))
    settransform!(vwalls[:left], Translation(0.0, d-wallThickness, 0.0))

    setobject!(vwalls[:right], MeshObject(
        Rect(Vec(0.0, 0.0, wallBottomEnd), Vec(0.0, wallThickness, wallTopEnd-wallBottomEnd)),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 1.0))))
    settransform!(vwalls[:right], Translation(0.0, d+D, 0.0))


    return vcart, vpendulum
end

function animate(Z)

    x1, x2 = Z[1][1:2]
    vcart, vpendulum = createAnimateObject(x1, x2)
    for z in Z[1:20:end]
        x1, x2 = z[1:2]
        settransform!(vcart, Translation(-0.2, x1-w/2.f0, 0.0))
        settransform!(vpendulum[:link], Translation(0.0, x1, 0.0) ∘ LinearMap(RotX(x2)))
        settransform!(vpendulum[:bob], Translation(0.0, x1, 0.0) ∘ LinearMap(RotX(x2)))
        sleep(0.04)
    end

end
