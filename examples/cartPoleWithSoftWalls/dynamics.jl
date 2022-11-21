using PyPlot 
using ControlSystems 

export CartPoleWithSoftWalls

const mc           = 1.0f0    
const mp           = 0.5f0 
const l            = 0.5f0        #wheel
const I1           = mp*l^2
const g            = 9.81f0
const d            = -0.5f0 
const D            = 1.0f0
const ϵn_const     = 1.0f0*ones(Float32, 2)
const ϵt_const     = 0.0f0*ones(Float32, 2)
const μ_const      = 0.0f0*ones(Float32, 2)
const gThreshold   = 0.001f0
const satu         = 100.0f0 

struct CartPoleWithSoftWalls{}  
end

function initialState(x1)
    # @assert pi-α <= θ0 <= pi+α "Give an initial spoke angle for the spoke in contact. This will help set the rimless wheel in contact with the surface"
    
    x1 = x1
    x2 = 5.0f0*pi/180f0

    x1dot = 0.0f0
    x2dot = 0.0f0

    return [x1, x2, x1dot, x2dot]
end

#returns the attributes needed to model contact
function (sys::CartPoleWithSoftWalls)(x, θ::Vector{T}; ϵn=ϵn_const, ϵt=ϵt_const, μ=μ_const, gThreshold=gThreshold, expert=false) where {T<:Real}
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
    g1 = q[1] - l*sin(q[2]) - d
    g2 = D - g1

    # return [g1, g2] 
    return [Inf, Inf]
end

function wn(z::Vector{T}) where {T<:Real}

    q, _ = parseStates(z)
    Wn = [1.0f0 -1.0f0;
        -l*cos(q[2]) l*cos(q[2])]

   return Wn 

end

function wt(z::Vector{T}) where {T<:Real}
    q, _ = parseStates(z)
    Wt = [0.0f0 0.0f0;
         -l*sin(q[2]) -l*sin(q[2])]
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
    Q = 10.0f0*Matrix{Float32}(LinearAlgebra.I, (4, 4)) 
    R = 3.0f0
    ControlSystems.lqr(cartLinearized, Q, R)
end

function lqr(z)
    k = vec([-1.82603  50.4554  -4.07912  15.4036])
    # k = zeros(4)
    return -k'*z
end

inputLayer(z) = [z[1], cos(z[2]), sin(z[2]), z[3], z[4]]

function control(z, θp::Vector{T}; expert=false) where {T<:Real}
    q, v = parseStates(z)

    if expert #working expert controller
        return clamp(lqr(z), -satu, satu)
    else
        return clamp(u, -satu, satu)
    end
end

function genForces(z, param::Vector{T}; expert=false) where {T<:Real}

    q, v = parseStates(z)
    x1, x2 = q
    #h = Bu - C qdot - G
    B = [1.0f0, 0.0f0]

    C = [0.0f0 mp*l*sin(x2); 
        -mp*sin(x2)/2.0f0 0.0f0]

    G = [0.0f0;
        -mp*g*l*sin(x2)]

    return B*control(z, param; expert=expert) - C*v - G
end

function createAnimateObject(x1, x2)
    vcart = vis[:cart]

    setobject!(vcart, MeshObject(
        Rect(Vec(0.0, 0.0, 0.0), Vec(0.2, 0.2, 0.1)),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 0.0, 0.25))))
    settransform!(vcart, Translation(-0.1, x1-0.1, 0.0))

    vpendulum = vis[:pendulum]
    setobject!(vpendulum[:link], MeshObject(
        Cylinder(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, l), 0.005),
        MeshLambertMaterial(color=RGBA{Float32}(1.0, 0.0, 0.0, 1.0))))
    settransform!(vpendulum[:link], Translation(0.0, x1, 0.0) ∘ LinearMap(RotX(x2)))

    setobject!(vpendulum[:bob], MeshObject(
        HyperSphere(Point(0.0, 0.0, l), 0.015),
        MeshLambertMaterial(color=RGBA{Float32}(0.0, 0.0, 1.0, 1.0))))
    settransform!(vpendulum[:bob], Translation(0.0, x1, 0.0) ∘ LinearMap(RotX(x2)))
   
    return vcart, vpendulum
end

function animate(Z)

    x1, x2 = Z[1][1:2]
    vcart, vpendulum = createAnimateObject(x1, x2)
    for z in Z[1:50:end]
        x1, x2 = z[1:2]
        settransform!(vcart, Translation(-0.1, x1-0.1, 0.0))
        settransform!(vpendulum[:link], Translation(0.0, x1, 0.0) ∘ LinearMap(RotX(x2)))
        settransform!(vpendulum[:bob], Translation(0.0, x1, 0.0) ∘ LinearMap(RotX(x2)))
        sleep(0.04)
    end

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
