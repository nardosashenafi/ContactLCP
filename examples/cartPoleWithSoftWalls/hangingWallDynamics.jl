using PyPlot 
using ControlSystems 
using MeshCat
using GeometryBasics
using CoordinateTransformations
using ColorTypes
using Blink
using Rotations

export CartPoleWithSoftWalls

##Posa's parameters
# const mc            = 1.0f0    
# const mp            = 0.5f0 
# const l             = 0.5f0     #center of mass of the pendulum  
# const I1            = mp*l^2.0f0/3.0f0
##Hardware parameters
const mc            = 0.57f0    
const mp            = 0.23f0 
const l             = 0.6413f0       #length of pendulum
const lcm           = l/2.0f0        #center of mass of the pendulum  
const I1            = mp*lcm^2.0f0/3.0f0
const g             = 9.81f0
const d             = -0.45f0 
const wallThickness = 0.10f0
const D             = 0.9f0
const wallBottomEnd = 0.45f0
const wallTopEnd    = 0.9f0
const contactNum    = 8
const ϵn_const      = 0.5f0*ones(Float32, contactNum)
const ϵt_const      = 0.0f0*ones(Float32, contactNum)
const μ_const       = 0.0f0*ones(Float32, contactNum)
const gThreshold    = 0.001f0
const satu          = 6.0f0     #Newtons. 10 Newton corresponds to 5.8 volts
const w             = 0.20f0

const leftWall_bottom   = [d, wallBottomEnd]
const leftWall_top      = [d, wallTopEnd]
const rightWall_bottom  = [D + d, wallBottomEnd]
const rightWall_top     = [D + d, wallTopEnd]
const lbr = leftWall_bottom     #bottom right corner of the left wall
const lbl = [leftWall_bottom[1] - wallThickness, leftWall_bottom[2]]  #bottom left corner of the left wall
const rbr = [rightWall_bottom[1] + wallThickness, rightWall_bottom[2]]  #bottom right corner of the right wall
const rbl = rightWall_bottom    #bottom left corner of the right wall

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
function (sys::CartPoleWithSoftWalls)(x, u::AbstractArray{T}; ϵn=ϵn_const, ϵt=ϵt_const, μ=μ_const, gThreshold=gThreshold, expert=false) where {T<:Real}
    q, v   = parseStates(x)
    x1,x2  = q[1:2]

    gn, Wn  = gap(x1, x2)  
    Wt  = wt(x1, x2)
    γn  = vnormal(Wn, v)
    γt  = vtang(Wt, v)
    M   = massMatrix(x2)
    h   = genForces(x, u; expert=expert)

    return gn, γn, γt, M, h, Wn, Wt, ϵn, ϵt, μ, gThreshold
end

#parses out the poses and velocities from the state vector
function (sys::CartPoleWithSoftWalls)(x::AbstractArray{T}) where {T<:Real}
    return parseStates(x)
end

function parseStates(x::AbstractArray{T}) where {T<:Real}
    q = x[1:2]      #[x1, x2]
    v = x[3:4]      #[x1dot, x2dot]
    return q, v
end

function pendulumPos(x, θ)
    return [x-l*sin(θ), l*cos(θ)]
end

function closestPointonPendulumToWall(pendulum_xy, x1, cornerx, cornery)
    px, py = pendulum_xy
    m = py/(px-x1)                     #pendulumSlope

    #draw circle around the corner and find the point on the pendulum (cx, cy) that is tangent to the circle
    cx = (m*cornery - cornerx + m^2*x1)/(m^2-1.0f0)
    cy = m*(cx-x1)

    return [cx, cy]
end

function gapPendulumToLeftWallCorners(pendulum_xy, x1, θ)

    c = closestPointonPendulumToWall(pendulum_xy, x1, lbr[1], lbr[2])
    g1inner = c[1] - lbr[1]

    ###also compute wn while you are here
    tt = tan(θ)
    A = (sec(θ))^2/(tt^2+1)^2*(lbr[2]*tt^2 - lbr[2] +2*lbr[1]*tt - 2*x1*tt)
    B = (sec(θ))^2/(tt^2+1)^2*(tt^2+1)

    wn1inner = [B; A]

    #now check the bottom left corners of the left wall
    c = closestPointonPendulumToWall(pendulum_xy, x1, lbl[1], lbl[2])
    g1outer = lbl[1] - c[1]
    A = (sec(θ))^2/(tt^2+1)^2*(lbl[2]*tt^2 - lbl[2] +2*lbr[1]*tt - 2*x1*tt)

    wn1outer = -[B; A]

    return g1inner, g1outer, wn1inner, wn1outer
end

function gapPendulumToRightWallCorners(pendulum_xy, x1, θ)

    c = closestPointonPendulumToWall(pendulum_xy, x1, rbl[1], rbl[2])
    g1inner = rbl[1] - c[1]

    ###also compute wn while you are here
    tt = tan(θ)
    A = (sec(θ))^2/(tt^2+1)^2*(lbr[2]*tt^2 - rbl[2] +2*rbl[1]*tt - 2*x1*tt)
    B = (sec(θ))^2/(tt^2+1)^2*(tt^2+1)

    wn1inner = -[B; A]

    #now check the bottom left corners of the left wall
    c = closestPointonPendulumToWall(pendulum_xy, x1, rbr[1], rbr[2])
    g1outer = c[1] - rbr[1]
    A = (sec(θ))^2/(tt^2+1)^2*(rbr[2]*tt^2 - rbr[2] +2*rbr[1]*tt - 2*x1*tt)

    wn1outer = [B; A]

    return g1inner, g1outer, wn1inner, wn1outer
end

function gapPendulumToLeftWall(pendulum_xy)
    g2innerleft = pendulum_xy[1] - d 
    g2outerleft = -g2innerleft - wallThickness
    wn2innerleft = [1.0f0; -l*cos(x2)]  
    wn2outerleft = -wn2innerleft

    return g2innerleft, g2outerleft, wn2innerleft, wn2outerleft
end

function gapPendulumToRightWall(pendulum_xy)
    g4innerleft = (D + d) - pendulum_xy[1]
    g4outerleft = -g4innerleft - wallThickness
    wn4innerleft = [-1.0f0; l*cos(x2)]   
    wn4outerleft = -wn4innerleft

    return g4innerleft, g4outerleft, wn4innerleft, wn4outerleft
end

function gap(x1, x2)
    #this function checks for two forms of contact on the two sides of the two walls (total of 8 cases)
    #The first set of gap functions checks if the pendulum link is in contact with the corners of the walls
    #The second set of gap functions check if the edge of the pendulum is contact with the sides of the walls
    pendulum_xy = pendulumPos(x1, x2)

    if pendulum_xy[2] > lbr[2] #pendulum in contact with left wall
        #check the length of the pendulum is in contact with the bottom corners of the wall
        g1innerleft, g1outerleft, wn1innerleft, wn1outerleft = gapPendulumToLeftWallCorners(pendulum_xy, x1, x2)                         
        #check if the end of the pendulum is in contact with the vertical wall sides
        g2innerleft, g2outerleft, wn2innerleft, wn2outerleft = gapPendulumToLeftWall(pendulum_xy)
    else
        #the pendulum is nowhere near the wall. Neither the length of the link nor its top end are in contact
        g1innerleft = norm(pendulum_xy - lbr, 2)
        g2innerleft = g1innerleft
        g1outerleft = -g1innerleft - wallThickness
        g2outerleft = g1outerleft
        x̄l, ȳl = pendulum_xy - lbr

        wn1innerleft = 1.0f0/g1innerleft*[x̄l; -x̄l*l*cos(x2) + ȳl*l*sin(x2)] 
        wn1outerleft = -wn1innerleft
        wn2innerleft = wn1innerleft
        wn2outerleft = -wn2innerleft
    end

    if pendulum_xy[2] > rbl[2] #pendulum in contact with right wall
        #check the length of the pendulum is in contact with the bottom corners of the wall
        g3innerleft, g3outerleft, wn3innerleft, wn3outerleft = gapPendulumToRightWallCorners(pendulum_xy, x1, x2)  
         #check if the end of the pendulum is in contact with the vertical wall sides
        g4innerleft, g4outerleft, wn4innerleft, wn4outerleft = gapPendulumToRightWall(pendulum_xy)
    else
        g3innerleft = norm(pendulum_xy - rbl, 2)
        g4innerleft = g3innerleft
        g3outerleft = -g3innerleft - wallThickness
        g4outerleft = g4innerleft

        x̄r, ȳr = pendulum_xy - rbl
        wn3innerleft = 1.0f0/g3innerleft*[x̄r; -x̄r*l*cos(x2) + ȳr*l*sin(x2)] 
        wn3outerleft = -wn3innerleft
        wn4innerleft = wn3innerleft
        wn4outerleft = -wn4innerleft
    end
    
    gn = [g1innerleft, g1outerleft, g2innerleft, g2outerleft, g3innerleft, g3outerleft, g4innerleft, g4outerleft] 
    wn = hcat(wn1innerleft, wn1outerleft, wn2innerleft, wn2outerleft, wn3innerleft, wn3outerleft, wn4innerleft, wn4outerleft)
    return gn, wn
end

function wt(x1, x2) 
    w1 = [0.0f0; -l*sin(x2)]

    return hcat(-w1, w1, -w1, w1, w1, -w1, w1, -w1)
end

function vnormal(Wn, v)
    return Wn'*v
end

function vtang(Wt, v)
    return Wt'*v
end

function massMatrix(x2)

    return [mc+mp -mp*lcm*cos(x2);
            -mp*lcm*cos(x2) mp*lcm^2+I1]
end

function lqrGains()
    M    = massMatrix(0.0f0)
    Minv = inv(M)

    Ĝ = Minv* [0.0f0;
              mp*g*lcm]
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

function lqr(z::AbstractArray{T}) where {T<:Real}
    k = convert.(T, vec([ -2.58192  32.2164  -4.41508  6.84527]))
    # k = zeros(4)
    return -k'*z
end

inputLayer(z) = [z[1], cos(z[2]), sin(z[2]), z[3], z[4]]

function control(z, u::AbstractArray{T}; expert=false, lqr_max = 10.0f0) where {T<:Real}
    q, v = parseStates(z)

    if expert #working expert controller
        return clamp(lqr(z), -lqr_max, lqr_max)
    else
        x1, x2 = q 
        x1dot, x2dot = v
        if ((1.0f0-cos(x2) <= (1.0f0-cosd(15.0))) && (abs(x2dot) <= 0.5f0))
            return clamp(lqr(z), -lqr_max, lqr_max)
        else
            return clamp(u[1], -satu, satu)
        end
    end
end

function genForces(x, u::AbstractArray{T}; expert=false) where {T<:Real}

    q, v = parseStates(x)
    x1, x2 = q[1:2]
    # B = [1.0f0, 0.0f0]

    # C = [0.0f0 mp*l*sin(x2); 
    #     -mp*sin(x2)/2.0f0 0.0f0]

    # G = [0.0f0;
    #     -mp*g*l*sin(x2)]

    #h = Bu - C qdot - G

    return  [1.0f0, 0.0f0]*control(x, u; expert=expert) -  #Bu 
                [0.0f0 mp*lcm*sin(x2); -mp*sin(x2)/2.0f0 0.0f0]*v -  #-C qdot
                [0.0f0; -mp*g*lcm*sin(x2)]            #-G
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
        Rect(Vec(0.0, 0.0, 0.0), Vec(0.2, w, 0.1)),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.25))))
    settransform!(vcart, Translation(-0.2, x1-w/2.0f0, 0.0))

    vpendulum = vis[:pendulum]
    setobject!(vpendulum[:link], MeshObject(
        Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, l), 0.005),
        MeshLambertMaterial(color=RGBA{Float32}(1.0, 0.0, 0.0, 1.0))))
    settransform!(vpendulum[:link], Translation(0.0, x1, 0.00) ∘ LinearMap(RotX(x2)))

    # setobject!(vpendulum[:bob], MeshObject(
    #     HyperSphere(Point(0.0, 0.0, l), 0.015),
    #     MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 1.0, 1.0))))
    # settransform!(vpendulum[:bob], Translation(0.0, x1, 0.00) ∘ LinearMap(RotX(x2)))
   
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
        # settransform!(vpendulum[:bob], Translation(0.0, x1, 0.0) ∘ LinearMap(RotX(x2)))
        sleep(0.04)
    end

end
