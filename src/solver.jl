
function lexiMin(z1, z2, LEXITHRESHOLD)   # picks the vector with the first zero or strictly negative index
    δz = z1 - z2

    if isempty(δz[abs.(δz) .> LEXITHRESHOLD])
        # println("vectors are equal! pick one")
        minz = (z1,1)
    else
        first(δz[abs.(δz) .> LEXITHRESHOLD]) < 0.0 ? minz = (z1,1) : minz = (z2,2)
    end

    return minz
end

function lexiMax(z1, z2)   # picks the vector with the first zero or strictly positive index
    δz = z1 - z2
    first(δz[abs.(δz) .> 1e-10]) < 0.0 ? maxz = (z2,2) : maxz = (z1,1)
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

function lexiblockRatioTest(Q::Matrix{T}, m; LEXITHRESHOLD=1e-10, zero_tol = 0.0, piv_tol = 1e-5 ) where {T<:Real}
    minIndex    = 1
    m_idIndex   = Vector{Int}() 
    #piv_tol  prevents it from selecting a pivot element that is almost zero

    [m[i] < -piv_tol ? push!(m_idIndex, i) : nothing for i in eachindex(m)]
    minIndex = first(m_idIndex)

    for i in 1:length(m_idIndex)-1
        j = m_idIndex[i+1]
        minlocalIndex = lexiMin((-Q[minIndex,:] .- zero_tol)./ m[minIndex], 
                                (-Q[j,:] .- zero_tol) ./ m[j], 
                                LEXITHRESHOLD)[2] 
        minlocalIndex == 2 ? minIndex = j : nothing 
    end 

    return minIndex  
end

function lexiPivot(Q, M, q::AbstractArray{T}, r, s) where {T<:Real}

    totalRow = range(1, stop=size(M, 1), step=1)
    totalCol = range(1, stop=size(M, 2), step=1)
    i = totalRow[totalRow .!= r]
    j = totalCol[totalCol .!= s]

    #copy current state to update with new states
    M̂ = T.(deepcopy(M))
    Q̂ = T.(deepcopy(Q))
    q̂ = T.(deepcopy(q))

    if abs.(M[r,s]) < 1e-5 
        println("About to divide by ", M[r,s])
    end

    q̂[r]    = -q[r]/M[r, s]
    q̂[i]    = q[i] - (M[i, s]/M[r,s])*q[r]
    M̂[r, s] = 1/M[r,s]
    M̂[i, s] = M[i, s]/M[r,s]
    M̂[r,j]  = -M[r,j]/M[r,s]
    M̂[i,j]  = M[i,j] - (M[i,s]/M[r,s])*M[r,j]'

    Q̂[r,:]  = -Q[r,:]/M[r,s]
    Q̂[i,:]  = Q[i,:] - (M[i,s]/M[r,s])*Q[r,:]'

    return Q̂, M̂, q̂
end

function lexiPivot!(Q, M, q::AbstractArray{T}, r, s) where {T<:Real}

    totalRow = range(1, stop=size(M, 1), step=1)
    totalCol = range(1, stop=size(M, 2), step=1)
    i = totalRow[totalRow .!= r]
    j = totalCol[totalCol .!= s]

    #copy current state to update with new states

    if abs.(M[r,s]) < 1e-5 
        println("About to divide by ", M[r,s])
    end

    #beware of the order of the following mutation
    q[i]    = q[i] - (M[i, s]/M[r,s])*q[r]
    q[r]    = -q[r]/M[r, s]
   
    Q[i,:]  = Q[i,:] - (M[i,s]/M[r,s])*Q[r,:]'
    Q[r,:]  = -Q[r,:]/M[r,s]
    
    M[i,j]  = M[i,j] - (M[i,s]/M[r,s])*M[r,j]'
    M[i, s] = M[i, s]/M[r,s]
    M[r,j]  = -M[r,j]/M[r,s]
    M[r, s] = 1/M[r,s]

end


function updateComplementPair(basic, nonbasic, oldBasicInd, newBasicIndex)
    locateComplementOfBasic = basic[oldBasicInd]
    nonbasic[locateComplementOfBasic] = newBasicIndex

    return nonbasic
end

function switchBasicWithNonBasic(basic, nonbasic, row, col)
    nonbasic = updateComplementPair(basic, nonbasic, row, col)
    basic[row], nonbasic[col] = nonbasic[col], basic[row]
    return basic, nonbasic
end

function complementPairIndex(nonbasic, nonbasicIndex)
    return nonbasic[nonbasicIndex]  #find complement of the dropped variable so we can add it to the basic variables on the next iteration
end

function lemkeLexi(M, q::AbstractArray{T}, x; MAX_ITER = 30, piv_tol = 1e-5) where {T<:Real}
    totalRow = size(M, 1)
    totalCol = size(M, 2)

    #starting values
    #z = -(M\q)

    if all(q .>= 0)
        return zeros(totalCol)
    end

    LEXITHRESHOLD   = 1e-10
    infeasibleCount = 0
    MAX_INFEASIBLECOUNT = 5
    iter            = 1

    pivottedIndices = Vector{Tuple{Int64, Int64}}()
    q̂               = deepcopy(q)

    while infeasibleCount < MAX_INFEASIBLECOUNT       #the value of LEXITHRESHOLD can cause infeasibility. In such cases, increase it.

        c  = ones(T, totalCol)    
        #z0 = maximum(-q ./ c)

        M̂ = [M c]

        Q = Matrix{T}(I, totalRow, totalRow)
        Q̂ = hcat(q̂, Q)

        basic    = collect(range(1, stop = totalRow, step= 1))
        nonbasic = collect(range(1, stop = totalCol+1, step= 1))

        alpha       = lexiαRatioTest(Q̂, c) 
        lexiPivot!(Q̂, M̂, q̂, alpha, totalCol+1)    #first feasible basic solution
        
        #keep track of the complementary pair indices
        basic, nonbasic = switchBasicWithNonBasic(basic, nonbasic, alpha, totalCol+1)

        d = alpha
        isFound         = false
        infeasible      = false
        iter            = 1

        while !isFound && !infeasible && iter < MAX_ITER

            if any(M̂[:,d] .< -piv_tol)
                b = lexiblockRatioTest(Q̂, M̂[:, d]; LEXITHRESHOLD=LEXITHRESHOLD)
            else
                LEXITHRESHOLD   *= 10   #when infeasible, make sure its not a thresholding issue
                infeasible      = true
                infeasibleCount += 1
                b = 1                   
                break           
            end

            #Pivotting
            if b == alpha
                lexiPivot!(Q̂ ,M̂, q̂, b, d)
                push!(pivottedIndices, (b, d))
                isFound = true
                infeasibleCount = 10
            else
                # This does not solve it at once
                lexiPivot!(Q̂ ,M̂, q̂, b, d)
                push!(pivottedIndices, (b, d))
                iter += 1
            end

            if !isFound
                basic, nonbasic = switchBasicWithNonBasic(basic, nonbasic, b, d)
                d = complementPairIndex(nonbasic, d) 
            end

        end     #end while
       
        if infeasibleCount < MAX_INFEASIBLECOUNT    #reset for another iteration
            pivottedIndices = Vector{Tuple{Int64, Int64}}()
            q̂               = deepcopy(q)
        end
    end     #end count to threshold adjustment
    
    if iter == MAX_ITER
        println("Exceeded max iteration. Increase your guess for q0")
    end

    if infeasibleCount == MAX_INFEASIBLECOUNT
        println("Infeasible")
        println("x = ", x)
    end

    # the solution is found in the vector q̂. 
    #But we have been switching the components of the vector.
    #so we need to trace back the exchange and extract the basic solution

    # basicSol = extractSolution(pivottedIndices, q̂, totalCol)
    basicSol = zeros(T, totalCol)
    return basicSol
end

function extractSolution(pivottedIndices, q̂::AbstractArray{T}, totalCol) where {T<:Real}

    basicSol = zeros(T, totalCol+1)
    pivotLen = length(pivottedIndices)
    windex  = getindex.(pivottedIndices, 1)
    zindex  = getindex.(pivottedIndices, 2)
    state_i = pivotLen

    for i in range(1, stop = totalCol+1, step=1 )
        extractionComplete = false
        basicSolRow = i
        m = findlast(x -> x == basicSolRow, windex[1:state_i])
        while !extractionComplete
            if isnothing(m)
                basicSol[i] = 0.0 
                extractionComplete = true
                break
            else
                state_i = m-1
                n = findlast(x -> x == zindex[m], zindex[1:state_i])
            end   

            if isnothing(n)
                basicSol[zindex[m]] = q̂[i] 
                extractionComplete = true
                break
            else
                basicSolRow = windex[n]
                # basicSol[zindex[n]] = q̂[windex[n]]    #q̂[windex[n]] belongs to w vector; not needed
                state_i = n-1 
            end
            m = findlast(x -> x == basicSolRow, windex[1:state_i])

            if isnothing(m)
                break 
            end
        end
        state_i = pivotLen
    end
    return basicSol[1:end-1]
end

#### the following functions are there for debugging purposes
function pivot(M, q::AbstractArray{T}, w, z, r, s) where {T<:Real}
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

function αRatioTest(q::AbstractArray{T}, d) where {T<:Real}
    ratio = Vector{T}(undef, length(q))
    [q[i] < 0.0 ? ratio[i] = -q[i]/d[i] : ratio[i] = Inf for i in eachindex(q)] 
    return argmin(ratio) + 1        # q and d are vectors with 3 components. In the augemented form, they are shifted by 1
end

function blockRatioTest(q::AbstractArray{T}, m) where {T<:Real}
    zero_tol = 1e-5
    ratio = Vector{T}(undef, length(m))
    [m[i]+zero_tol < 0.0 ? ratio[i] = -q[i]/m[i] : ratio[i] = Inf for i in eachindex(m)] 
    return argmin(ratio) 
end

function lemke(M, q::AbstractArray{T}) where {T<:Real}
    # println("Lemke begins")
    totalRow = size(M, 1)
    totalCol = size(M, 2)

    #starting values
    z = -(M\q)
    w = zeros(totalRow)

    if all(q .>= 0)
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

function stepLemke(M, q::AbstractArray{T}, x) where {T<:Real}
    # println("x = ", x)

    totalRow = size(M, 1)
    totalCol = size(M, 2)

    #starting values
    z = -(M\q)
    w = zeros(totalRow)

    if all(q .>= 0)
        return zeros(totalCol)
    end

    #Trivial solution does not apply. Create the augemented form
    pivottedIndices = Vector{Tuple{Int64, Int64}}()
    c  = ones(totalCol)    
    q0 = 10000.0f0              #start with sufficiently large scalar q0 >= 0.0
    w0 = 0.0
    z0 = maximum(-q ./ c)

    ẑ = vcat(z0, z)
    # ŵ = vcat(w0, w)
    # q̂ = vcat(q0, q)
    ŵ = deepcopy(w)
    q̂ = deepcopy(q)
    M̂ = [M c]

    Q = Matrix{Float32}(I, totalRow, totalRow)
    # Q̂ = zeros(totalRow+1, totalCol+1)
    # Q̂[:,1] = q̂
    # Q̂[2:end, 2:end] = Q     # Q̂ = [q0 0; q Q]
    Q̂ = hcat(q̂, Q)

    # basic = [1, 2, 3]   #the complement pair indices of the basic variables in the nonbasic variables
    # nonbasic = [1, 2, 3, -1]    #the complement pair indices of the nonbasic variables in the basic variables

    basic = collect(range(1, stop = totalRow, step= 1))
    nonbasic = collect(range(1, stop = totalCol+1, step= 1))
    # basicDict = Dict("w1" => 1, "w2" => 2, "w3" => 3)
    # nonbasicDict = Dict("z1" => 1, "z2" => 2, "z3" => 3, "z0" => -1)

    α = lexiαRatioTest(Q̂, c) 
    Q̂, M̂, q̂, ŵ, ẑ = lexiPivot(Q̂, M̂, q̂, ŵ, ẑ, α, totalCol+1)    #first feasible basic solution
    basic, nonbasic = switchBasicWithNonBasic(basic, nonbasic, α, totalCol+1)
    # basicDict, nonbasicDict = switchComplementInDict(basicDict, nonbasicDict, α, totalCol+1)
    # push!(pivottedIndices, (α, 1))
    # println("α = ", α)
    # Minitial = deepcopy(M̂)
    # Minitial[:, 1] = zeros(totalRow)
    # Minitial[α, 1] = 1.0

    B = Matrix{Float32}(I, totalRow, totalRow)
    # Binv = hcat(q̂, inv(B))
    d = α
    isFound     = false
    infeasible  = false
    MAX_ITER    = 30
    iter        = 1

    piv_tol = 1e-5  
    while !isFound && !infeasible && iter < MAX_ITER

        # println("m = ", M̂[:,d])
        # println("q̂ = ", q̂)
        #determine blocking variable
        if any(M̂[:,d] .< -piv_tol)
            # m = deepcopy((M̂[:,d]))
            # minIndex = 1
            # m_idIndex = Vector{Int}() 
            # piv_tol = 0.0
            # zero_tol = 0.0
            # [m[i]+piv_tol < 0.0 ? push!(m_idIndex, i) : nothing for i in 1:length(m)]
        
            # println("m = ", m)
            # println("Relevant Q ratio ", [(-Q̂[j,:] .- zero_tol)./ m[j] for j in m_idIndex])

            # println("Binv = ", Binv)
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
        # println("index = ", pivottedIndices[end])

        if !isFound
            basic, nonbasic = switchBasicWithNonBasic(basic, nonbasic, b, d)
            d = complementPairIndex(nonbasic, d) 
        end
        # basicDict, nonbasicDict = switchComplementInDict(basicDict, nonbasicDict, b, d)

              # driving variable is the complement of the current blocking variable
        # r, s = pivottedIndices[end]
        # B[:,r] = -Minitial[:,s] 
        # Binv = hcat(q̂, inv(B))
    end     #end while
    
    if iter == MAX_ITER
        println("Exceeded max iteration. Increase your guess for q0")
    end


    basicSol = zeros(T, totalCol+1)
    pivotLen = length(pivottedIndices)
    windex  = getindex.(pivottedIndices, 1)
    zindex  = getindex.(pivottedIndices, 2)
    state_i = pivotLen

    for i in range(1, stop = totalCol+1, step=1 )
        extractionComplete = false
        m = findlast(x -> x == i, windex[1:state_i])
        while !extractionComplete
            if isnothing(m)
                basicSol[i] = 0.0 
                extractionComplete = true
                break
            else
                state_i = m-1
                n = findlast(x -> x == zindex[m], zindex[1:state_i])
            end   

            if isnothing(n)
                basicSol[zindex[m]] = q̂[windex[m]] 
                extractionComplete = true
                break
            else
                i = windex[n]
                # basicSol[zindex[n]] = q̂[windex[n]]    #q̂[windex[n]] belongs to w vector; not needed
                state_i = n-1 
            end
            m = findlast(x -> x == i, windex[1:state_i])

            if isnothing(m)
                break 
            end
        end
        state_i = pivotLen
    end

    basicSol = basicSol[1:end-1]
    # println("Pivotted indices = ", pivottedIndices)

    solpathsolver = lcpOpt(M, q, Int(floor(length(q)/3)))

    if any(abs.(basicSol - solpathsolver) .> 0.001)
        println("Pivotted indices = ", pivottedIndices)

        println("pathsolver = ", solpathsolver)
        println("Lemke = ", basicSol)
        println("q̂ = ", q̂)
        println("x = ", x)
    end
    return basicSol
end