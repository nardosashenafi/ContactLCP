
function pivot(M, q, w, z, r, s)
    totalRow = range(1, stop=length(w), step=1)
    totalCol = range(1, stop=length(z), step=1)
    i = totalRow[totalRow .!= r]
    j = totalCol[totalCol .!= s]

    #copy current state to update with new states
    ŵ = deepcopy(w)
    ẑ = deepcopy(z)
    M̂ = deepcopy(M)
    q̂ = deepcopy(q)

    ŵ[r] = z[s]
    ẑ[s] = w[r]
    q̂[r] = -q[r]/M[r, s]
    q̂[i] = q[i] - (M[i, s]/M[r,s])*q[r]
    M̂[r, s] = 1/M[r,s]
    M̂[i, s] = M[i, s]/M[r,s]
    M̂[r,j] = -M[r,j]/M[r,s]
    M̂[i,j] = M[i,j] - (M[i,s]/M[r,s])*M[r,j]'

    return M̂, q̂, ŵ, ẑ
end

function αRatioTest(q, d)
    ratio = Vector{Float32}(undef, length(q))
    [q[i] < 0.0 ? ratio[i] = -q[i]/d[i] : ratio[i] = Inf for i in eachindex(q)] 
    return argmin(ratio) + 1        # q and d are vectors with 3 components. In the augemented form, they are shifted by 1
end

function blockRatioTest(q, m)
    zero_tol = 1e-5
    ratio = Vector{Float32}(undef, length(m))
    [m[i]+zero_tol < 0.0 ? ratio[i] = -q[i]/m[i] : ratio[i] = Inf for i in eachindex(m)] 
    return argmin(ratio) 
end

function lemke(M, q)
    println("Lemke begins")
    totalRow = size(M, 1)
    totalCol = size(M, 2)

    #starting values
    z = -(M\q)
    w = [0.0, 0.0, 0.0]

    if all(q .>= 0)
        println("trivial solution")
        w = q 
        println("w = ", w)
        return w
    end

    #Trivial solution does not apply. Create the augemented form
    pivottedIndices = Vector{Tuple{Int64, Int64}}()
    q0 = 2.0f0                #start with sufficiently large scalar q0 >= 0.0
    c  = 1.0f0*ones(totalCol)                  # c > 0
    # w0 = q0 - c'*z                 # solves original LCP when z0 = 0
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
    MAX_ITER    = 10
    iter        = 1

    while !isFound && !infeasible && iter < MAX_ITER
        #determine blocking variable
        # println("M̂ = ", M̂)
        if any(M̂[:,d] .< 0.0f0)
            b = blockRatioTest(q̂, M̂[:, d])
            # println("blocks = ", b, " driving = ", d)
        else
            println("Interpret output interms of infeasibility or unsolvability")
            infeasible = true
            b = 1                   #w0 is the blocking variable
            break            #TODO: handle infeasibility
        end

        #Pivotting
        if b == α
            M̂, q̂, ŵ, ẑ = pivot(M̂, q̂, ŵ, ẑ, b, d)
            push!(pivottedIndices, (b, d))
            println("Solved")
            # println("q̂ = ", q̂)
            isFound = true

        else
            # This does not solve it at once
            M̂, q̂, ŵ, ẑ = pivot(M̂, q̂, ŵ, ẑ, b, d)    # z_d = z0
            push!(pivottedIndices, (b, d))
            # println("q̂ = ", q̂)
            iter += 1
            # println("iterating")
        end
        d = b
    end     #end while
    

    if iter == MAX_ITER
        println("Exceeded max iteration. Increase your guess for q0")
    end

    println(pivottedIndices)

    sol = zeros(totalCol)
    for (r, s) in pivottedIndices
        sol[s-1] = q̂[r]
    end

    # println("q̂ = ", q̂)
    println("pathsolver = ", lcpOpt(M, q, 1))
    println("Lemke = ", sol)
    return sol
end
