
function parseStates(x) 
    q = x[1:2]
    u = x[3:4]
    return q, u
end

function impactMap(x)

    ϕ = x[2]
    det = I1*I2 + I1*m2*l2^2 + I2*mt*l1^2 + 
            m2*l1^2.0f0*l2^2.0f0*(m1 + m2*(sin(α - ϕ))^2.0f0)

    ξ1 = 1.0f0/det * ((I1*I2 + I1*m2*l2^2.0f0) + 
            (I2*mt*l1^2.0f0 + m2*l1^2.0f0*l2^2.0f0*(m1 +0.5f0* m2))*cos(2.0f0α) -
            1.0f0/2.0f0*m2^2.0f0*l1^2.0f0*l2^2.0f0*cos(2.0f0ϕ))

    ξ2 = 1.0f0/det *(m2*l1*l2*(I1*(cos(α - ϕ) - cos(α + ϕ)) + 
            mt*l1^2.0f0*(cos(2.0f0α)*cos(α - ϕ) - cos(α + ϕ))) )

    return [ξ1 0.0f0;
            ξ2 1.0f0]

end

function impactMap!(Ξ, x)

    ϕ = x[2]
    det = I1*I2 + I1*m2*l2^2 + I2*mt*l1^2 + 
            m2*l1^2.0f0*l2^2.0f0*(m1 + m2*(sin(α - ϕ))^2.0f0)

    ξ1 = 1.0f0/det * ((I1*I2 + I1*m2*l2^2.0f0) + 
            (I2*mt*l1^2.0f0 + m2*l1^2.0f0*l2^2.0f0*(m1 +0.5f0* m2))*cos(2.0f0α) -
            1.0f0/2.0f0*m2^2.0f0*l1^2.0f0*l2^2.0f0*cos(2.0f0ϕ))

    ξ2 = 1.0f0/det *(m2*l1*l2*(I1*(cos(α - ϕ) - cos(α + ϕ)) + 
            mt*l1^2.0f0*(cos(2.0f0α)*cos(α - ϕ) - cos(α + ϕ))) )

    Ξ[1,1] = ξ1
    Ξ[2,1] = ξ2

end

function gap(x)
    return [l1*cos(x[1]) - l1*cos(2.0f0α - abs(x[1]))] #l1cos(θ) - l1*cos(2α - |θ|)
end

function gap!(gn, x)
    gn[1] = l1*cos(x[1]) - l1*cos(2.0f0α - abs(x[1]))
end

function massMatrix(x)
    θ, ϕ = x[1:2]
    M    = [I1 + mt*l1^2 -m2*l1*l2*cos(θ - ϕ);
            -m2*l1*l2*cos(θ - ϕ) I2 + m2*l2^2]
end

function massMatrix!(M, x)
    θ, ϕ = x[1:2]
    M[:,:]  = [I1 + mt*l1^2 -m2*l1*l2*cos(θ - ϕ);
                -m2*l1*l2*cos(θ - ϕ) I2 + m2*l2^2]
end

function vnormal(x)
    sgn = sign(x[1])
    _, u = parseStates(x)
    return (-l1*sin(x[1]) - sgn*l1*sin(2.f0α - abs(x[1])))*u[1] + 0.0f0*u[2]
end

function vnormal!(vn, x)
    sgn = sign(x[1])
    _, u = parseStates(x)
    vn[1] = (-l1*sin(x[1]) - sgn*l1*sin(2.f0α - abs(x[1])))*u[1] + 0.0f0*u[2]
end

wrap(x) = [cos(x[1]), sin(x[1]), cos(x[2]), sin(x[2]), x[3], x[4]]

function control(x, θp)
    @assert length(θp) == 6 + DiffEqFlux.paramlength(Hd)
    y = MLBasedESC.controller(npbc, wrap(x), θp)
    return clamp(y, -satu, satu)
end

function control!(τ, x, θp)
    y = MLBasedESC.controller(npbc, CuArray(wrap(x)), θp)
    τ[1] = clamp(y, -satu, satu)
end

function genForces(x)
    q, u = parseStates(x)
    θ, ϕ = q[1:2]
    #h = Bu - C qdot - G
    B = [-1.0f0, 1.0f0]
    C = m2*l1*l2*sin(θ - ϕ)*[0.0f0 -u[2];
                            u[1] 0.0f0]

    G = g*[-mt*l1*sin(θ - γ),
        m2*l2*sin(ϕ - γ)]

    return B, C, G
end

function genForces!(B, C, G, x)
    q, u = parseStates(x)
    θ, ϕ = q[1:2]
    #h = Bu - C qdot - G
    B[:] = [-1.0f0, 1.0f0] 
    C[:,:] = m2*l1*l2*sin(θ - ϕ)*[0.0f0 -u[2];
                                u[1] 0.0f0]

    G[:,:] = g*[-mt*l1*sin(θ - γ),
                m2*l2*sin(ϕ - γ)]

end

function checkContact(gn) 
     
    gn[i] < gThreshold ? contactIndex = [0.0f0] : contactIndex = [1.0f0] 
    return contactIndex
end

function checkContact!(contactIndex, gn) 
     
    for i in 1:total_contact_num
        if gn[i] < gThreshold 
            contactIndex[i] = 1.0f0 
        end
    end
end

function allocateCuArrays(x1, θ)
    gn          = CuArray(gap(x1))
    M           = CuArray(massMatrix(x1))
    B, C, G     = CuArray.(genForces(x1))

    Ξ           = CuArray(impactMap(x1))
    contactIndex  = CuArray(checkContact(gn))
    τ             = CUDA.zeros(1)
    control!(τ, x, θ)
    Y               = Vector{CuArray{Float32, 1}}()
    tY              = Vector{CuArray{Float32, 1}}()

    return gn, M, B, C, G, τ, Ξ, contactIndex, Y, tY
end

Zygote.@adjoint CUDA.zeros(x...) = CUDA.zeros(x...), _ -> map(_ -> nothing, x)

function fulltimestep(x1, θ; timeSteps = totalTimeStep) 
    
    Y  = Vector{Vector{Float32}}()
    tY = Vector{Float32}()
    y  = [deepcopy(x1)]
    ty = 0.0f0
    
    for i in 1:timeSteps
        qA, uA      = parseStates(x1)
        qM          = qA + 0.5f0*Δt*uA
        
        x_mid       = vcat(qM,uA)

        gn           = gap(x_mid)
        M            = massMatrix(x_mid)
        contactIndex = checkContact(gn)
        B, C, G = genForces(x_mid)
        h       = B*control(x_mid, θ) - C*uA - G

        if any(contactIndex .== 1.0f0) && any(vnormal(x_mid) .< 0.0f0)

            Ξ   = impactMap(x_mid)
            uE  = Ξ*uA
            #switch θ to the next spoke
            if uA[1] < 0.0
                q1M = α
            elseif uA[1] > 0.0
                q1M = -α
            else
                q1M = qM[1]
            end
            
            qE = [q1M, qM[2]] + 0.5f0*Δt*uE
            x2 = vcat(qE,uE)

            y   = vcat(y, [x2])
            ty  = vcat(ty, ty[end]+Δt)
            Y   = vcat(Y, [y])
            tY  = vcat(tY, [ty])

            y   = [deepcopy(y[end])]
            ty  = ty[end]
        else

            uE  = M\(h*Δt) + uA
            qE  = qM + 0.5f0*Δt*uE
            x2  = vcat(qE,uE)
            y   = vcat(y, [x2])
            ty  = vcat(ty, ty[end]+Δt)
        end
        x1 = x2
    end
    return Y, tY
end


function fulltimestep!(x1, param, gn, M, B, C, G, τ, Ξ, contactIndex, Y, tY; timeSteps = totalTimeStep) 
    
    y  = [deepcopy(x1)]
    ty = CuArray([0.0f0])
    
    for i in 1:timeSteps
        qA, uA      = parseStates(x1)
        qM          = qA + 0.5f0*Δt*uA
        
        x_mid       = vcat(qM,uA)

        gap!(gn, x_mid)
        massMatrix!(M, x_mid)
        checkContact!(contactIndex, gn)
        genForces!(B, C, G, x_mid)
        control!(τ, x_mid, param)
        h   = B*τ[1] - C*uA - G

        if any(contactIndex .== 1.0f0) && any(vnormal(x_mid) .< 0.0f0)

            impactMap!(Ξ, x_mid)
            uE = Ξ*uA
            #switch θ to the next spoke
            if uA[1] < 0.0
                qM[1] = α
            elseif uA[1] > 0.0
                qM[1] = -α
            end
            
            qE = qM + 0.5f0*Δt*uE
            x2 = vcat(qE,uE)

            y   = vcat(y, [x2])
            ty  = vcat(ty, ty[end]+Δt)
            Y   = vcat(Y, [y])
            tY  = vcat(tY, [ty])

            y   = [deepcopy(y[end])]
            ty  = ty[end]
        else

            uE  = M\(h*Δt) + uA
            qE  = qM + 0.5f0*Δt*uE
            x2  = vcat(qE,uE)
            y   = vcat(y, [x2])
            ty  = vcat(ty, ty[end]+Δt)
        end
        x1 = x2
    end
    Y   = vcat(Y, [y])
    tY  = vcat(tY, [ty])

    ##compute loss
    loss = hipSpeedLoss(Y, tY)
end

