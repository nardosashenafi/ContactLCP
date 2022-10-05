
function pivot(M, q::Vector{T}, w, z, r, s) where {T<:Real}
    totalRow = range(1, stop=length(w), step=1)
    totalCol = range(1, stop=length(z), step=1)
    i = totalRow[totalRow .!= r]
    j = totalCol[totalCol .!= s]

    #copy current state to update with new states
    ŵ = T.(deepcopy(w))
    ẑ = T.(deepcopy(z))
    M̂ = T.(deepcopy(M))
    q̂ = T.(deepcopy(q))

    ŵ[r]    = z[s]
    ẑ[s]    = w[r]
    q̂[r]    = -q[r]/M[r, s]
    q̂[i]    = q[i] - (M[i, s]/M[r,s])*q[r]
    M̂[r, s] = 1/M[r,s]
    M̂[i, s] = M[i, s]/M[r,s]
    M̂[r,j]  = -M[r,j]/M[r,s]
    M̂[i,j]  = M[i,j] - (M[i,s]/M[r,s])*M[r,j]'

    return M̂, q̂, ŵ, ẑ
end

function αRatioTest(q::Vector{T}, d) where {T<:Real}
    ratio = Vector{T}(undef, length(q))
    [q[i] < 0.0 ? ratio[i] = -q[i]/d[i] : ratio[i] = Inf for i in eachindex(q)] 
    return argmin(ratio) + 1        # q and d are vectors with 3 components. In the augemented form, they are shifted by 1
end

function blockRatioTest(q::Vector{T}, m) where {T<:Real}
    zero_tol = 1e-5
    ratio = Vector{T}(undef, length(m))
    [m[i]+zero_tol < 0.0 ? ratio[i] = -q[i]/m[i] : ratio[i] = Inf for i in eachindex(m)] 
    return argmin(ratio) 
end

function lemke(M, q::Vector{T}) where {T<:Real}
    # println("Lemke begins")
    totalRow = size(M, 1)
    totalCol = size(M, 2)

    #starting values
    z = -(M\q)
    w = zeros(totalRow)

    if all(q .>= 0)
        w = q 
        return zeros(totalCol)
    end

    #Trivial solution does not apply. Create the augemented form
    pivottedIndices = Vector{Tuple{Int64, Int64}}()
    c  = 1.0f0*ones(totalCol)                  # c > 0

    aug_size = 1
    q0 = 10000.0f0              #start with sufficiently large scalar q0 >= 0.0
    w0 = 0.0
    z0 = maximum(-q ./ c)

    ẑ = vcat(z0, z)
    ŵ = vcat(w0, w)
    q̂ = vcat(q0, q)
    M̂ = [0.0f0 -c'; c M]

    α = αRatioTest(q, c) 
    M̂, q̂, ŵ, ẑ = pivot(M̂, q̂, ŵ, ẑ, α, 1)    #first feasible basic solution
    # push!(pivottedIndices, (α, 1))

    d = α
    isFound     = false
    infeasible  = false
    MAX_ITER    = 100
    iter        = 1

    while !isFound && !infeasible && iter < MAX_ITER

        #determine blocking variable
        if any(M̂[:,d] .< 0.0f0)
            b = blockRatioTest(q̂, M̂[:, d])
        else
            println("Interpret output interms of infeasibility or unsolvability")
            infeasible = true
            b = 1                   
            break            #TODO: handle infeasibility
        end

        #Pivotting
        if b == α

            M̂, q̂, ŵ, ẑ = pivot(M̂, q̂, ŵ, ẑ, b, d)
            push!(pivottedIndices, (b, d))
            # println("Solved")
            isFound = true

        else
            # This does not solve it at once
            M̂, q̂, ŵ, ẑ = pivot(M̂, q̂, ŵ, ẑ, b, d)    # z_d = z0
            push!(pivottedIndices, (b, d))
            iter += 1
        end
        d = b               # driving variable is the complement of the current blocking variable
    end     #end while
    
    println("Pivotted indices = ", pivottedIndices)
    if iter == MAX_ITER
        println("Exceeded max iteration. Increase your guess for q0")
    end

    sol = zeros(T, totalCol)
    for (r, s) in pivottedIndices
        sol[s-1] = q̂[r]
    end

    # println("q̂ = ", q̂)
    println("pathsolver = ", lcpOpt(M, q, Int(floor(length(q)/3))))
    println("Lemke = ", sol)
    return sol
end

function lexiMin(z1, z2)   # picks the vector with the first zero or strictly negative index
    δz = z1 - z2
    first(δz[δz .!= 0.0]) < 0.0 ? minz = (z1,1) : minz = (z2,2)
    return minz
end

function lexiMax(z1, z2)   # picks the vector with the first zero or strictly positive index
    δz = z1 - z2
    first(δz[δz .!= 0.0]) < 0.0 ? maxz = (z2,2) : maxz = (z1,1)
    return maxz
end

function lexiαRatioTest(Q::Matrix{T}, d) where {T<:Real}

    Qi0Index = Vector{Int}() 
    [Q[i,1] < 0.0 ? push!(Qi0Index, i) : nothing for i in 1:size(Q, 1)]

    minIndex = first(Qi0Index)
    for i in 1:length(Qi0Index)-1
        j = Qi0Index[i+1]
        minlocalIndex = lexiMax(-Q[minIndex,:] ./ d[minIndex], -Q[j,:] ./ d[j])[2]
        minlocalIndex == 2 ? minIndex = j : nothing  
    end 

    return minIndex         # q and d are vectors with 3 components. In the augemented form, they are shifted by 1
end

function lexiblockRatioTest(Q::Matrix{T}, m) where {T<:Real}
    minIndex = 1
    m_idIndex = Vector{Int}() 
    [m[i] < 0.0 ? push!(m_idIndex, i) : nothing for i in 1:length(m)]

    minIndex = first(m_idIndex)
    for i in 1:length(m_idIndex)-1
        j = m_idIndex[i+1]
        minlocalIndex = lexiMin(-Q[minIndex,:] ./ m[minIndex], -Q[j,:] ./ m[j])[2] 
        minlocalIndex == 2 ? minIndex = j : nothing 
    end 

    return minIndex  
end

function lexiPivot(Q, M, q::Vector{T}, w, z, r, s) where {T<:Real}

    totalRow = range(1, stop=length(w), step=1)
    totalCol = range(1, stop=length(z), step=1)
    i = totalRow[totalRow .!= r]
    j = totalCol[totalCol .!= s]

    #copy current state to update with new states
    ŵ = T.(deepcopy(w))
    ẑ = T.(deepcopy(z))
    M̂ = T.(deepcopy(M))
    Q̂ = T.(deepcopy(Q))
    q̂ = T.(deepcopy(q))

    ŵ[r]    = z[s]
    ẑ[s]    = w[r]
    q̂[r]    = -q[r]/M[r, s]
    q̂[i]    = q[i] - (M[i, s]/M[r,s])*q[r]
    M̂[r, s] = 1/M[r,s]
    M̂[i, s] = M[i, s]/M[r,s]
    M̂[r,j]  = -M[r,j]/M[r,s]
    M̂[i,j]  = M[i,j] - (M[i,s]/M[r,s])*M[r,j]'

    Q̂[r,:]  = -Q[r,:]/M[r,s]
    Q̂[i,:]  = Q[i,:] - (M[i,s]/M[r,s])*Q[r,:]'

    return Q̂, M̂, q̂, ŵ, ẑ
end

function lemkeLexi(M, q::Vector{T}) where {T<:Real}
    # println("Lemke begins")
    totalRow = size(M, 1)
    totalCol = size(M, 2)

    # sol = zeros(T, totalCol+1)
    #starting values
    z = -(M\q)
    w = zeros(totalRow)

    if all(q .>= 0)
        println("Trivial solution = ", w)
        println("pathsolver = ", lcpOpt(M, q, Int(floor(length(q)/3))))
        return zeros(totalCol)
    end

    #Trivial solution does not apply. Create the augemented form
    pivottedIndices = Vector{Tuple{Int64, Int64}}()
    c  = 1.0f0*ones(totalCol)                  # c > 0
    q0 = 1000.0f0              #start with sufficiently large scalar q0 >= 0.0
    w0 = 0.0
    z0 = maximum(-q ./ c)

    ẑ = vcat(z0, z)
    ŵ = vcat(w0, w)
    q̂ = vcat(q0, q)
    M̂ = [0.0f0 -c'; c M]

    Q = Matrix{Float32}(I, totalRow, totalRow)
    Q̂ = zeros(totalRow+1, totalCol+1)
    Q̂[:,1] = q̂
    Q̂[2:end, 2:end] = Q     # Q̂ = [q0 0; q Q]
    # Q̂ = hcat(q̂, Q)

    α = lexiαRatioTest(hcat(q, Q), c) + 1
    Q̂, M̂, q̂, ŵ, ẑ = lexiPivot(Q̂, M̂, q̂, ŵ, ẑ, α, 1)    #first feasible basic solution
    # push!(pivottedIndices, (α, 1))
    # sol[α] = q̂[1]

    d = α
    isFound     = false
    infeasible  = false
    MAX_ITER    = 100
    iter        = 1

    while !isFound && !infeasible && iter < MAX_ITER

        #determine blocking variable
        if any(M̂[:,d] .< 0.0f0)
            b = lexiblockRatioTest(Q̂, M̂[:, d])
        else
            println("Interpret output interms of infeasibility or unsolvability")
            infeasible = true
            b = 1                   
            break            #TODO: handle infeasibility
        end

        #Pivotting
        if b == α

            Q̂, M̂, q̂, ŵ, ẑ = lexiPivot(Q̂ ,M̂, q̂, ŵ, ẑ, b, d)
            push!(pivottedIndices, (b, d))
            # println("Solved")
            isFound = true

        else
            # This does not solve it at once
            Q̂, M̂, q̂, ŵ, ẑ = lexiPivot(Q̂ ,M̂, q̂, ŵ, ẑ, b, d)
            push!(pivottedIndices, (b, d))
            iter += 1
        end
        # sol[b] = q̂[d]
        d = b               # driving variable is the complement of the current blocking variable
    end     #end while
    
    println("Pivotted indices = ", pivottedIndices)
    if iter == MAX_ITER
        println("Exceeded max iteration. Increase your guess for q0")
    end

    sol = zeros(T, totalCol)
    count_unique = []
    for (r, s) in reverse(pivottedIndices)
        r ∉ count_unique ? sol[s-1] = q̂[r] : nothing
        push!(count_unique, r)
    end
    # if totalCol > 4
    #     sol[3] = 0.0
    # end

    println("pathsolver = ", lcpOpt(M, q, Int(floor(length(q)/3))))
    println("Lemke = ", sol)
    return sol
end

