using PyPlot 
using ControlSystems 
using MeshCat
using GeometryBasics
using CoordinateTransformations
using ColorTypes
using Blink
using Rotations

export CartPoleWithSoftWalls

const mc           = 1.0f0    
const mp           = 0.5f0 
const l            = 0.5f0        #wheel
const I1           = mp*l^2
const g            = 9.81f0
const d            = -0.75f0 
const D            = 1.5f0
const ϵn_const     = 0.5f0*ones(Float32, 4)
const ϵt_const     = 0.0f0*ones(Float32, 4)
const μ_const      = 0.0f0*ones(Float32, 4)
const gThreshold   = 0.001f0
const satu         = 4.0f0
const w            = 0.20f0

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
    q = x[1:2]      #[x1, y2]
    u = x[3:4]      #[x1dot, x2dot]
    return q, u
end

function gap(z)
    q, _ = parseStates(z)    
    g1 = q[1] - l*sin(q[2]) - d #pendulum in contact with wall
    g2 = D - g1
    g3 = q[1] - w/2.0f0 - d     #cart in contact with wall
    g4 = D - g3 - w
    
    return [g1, g2, g3, g4] 
    # return [g1, g2] 
    # return [Inf, Inf, Inf, Inf]
end

function wn(z::Vector{T}) where {T<:Real}

    q, _ = parseStates(z)
    Wn = [1.0f0 -1.0f0 1.0f0 -1.0f0;
        -l*cos(q[2]) l*cos(q[2]) 0.0f0 0.0f0]
    # Wn = [1.0f0 -1.0f0;
    #     -l*cos(q[2]) l*cos(q[2])]

   return Wn 

end

function wt(z::Vector{T}) where {T<:Real}
    q, _ = parseStates(z)
    Wt = [0.0f0 0.0f0 0.0f0 0.0f0;
         -l*sin(q[2]) -l*sin(q[2]) 0.0f0 0.0f0]
    # Wt = [0.0f0 0.0f0;
    #     -l*sin(q[2]) -l*sin(q[2])]

    return Wt
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
        if ((1.0f0-cos(x2) <= 1.0f0-cosd(17.0)) && x2dot <= 0.5f0)
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
        Rect(Vec(0.0, 0.0, -0.8), Vec(0.0, 0.02, 1.6)),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 1.0))))
    settransform!(vwalls[:left], Translation(0.0, d, 0.0))

    setobject!(vwalls[:right], MeshObject(
        Rect(Vec(0.0, 0.0, -0.8), Vec(0.0, 0.02, 1.6)),
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
