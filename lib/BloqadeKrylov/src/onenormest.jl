
function onenormest_explicit(A, p::Int=1)
    if p <=0
       error("p must be positive") 
    end

    if p == 1
        return opnorm(A,1)
    else
        D = similar(A)
        mul!(D,A,A)
        for i in 3:p
            mul!(D,A,D)
        end
    end

    return opnorm(D,1)
end

"""
    Compute a lower bound of the 1-norm of a matrix to a power of p

    Parameters
    ----------
    A : ndarray or other linear operator
        A linear operator that can be transposed and that can
        produce matrix products.
    p : power of A, optional (default =1)
    t : int, optional
        A positive parameter controlling the tradeoff between
        accuracy versus time and memory usage.
        Larger values take longer and use more memory
        but give more accurate output.
    itmax : int, optional
        Use at most this many iterations.

    Returns
    -------
    est : float
        An underestimate of the 1-norm of the sparse matrix.

    Notes
    -----
    This is algorithm 2.4 of [1].

    In [2] it is described as follows.
    "This algorithm typically requires the evaluation of
    about 4t matrix-vector products and almost invariably
    produces a norm estimate (which is, in fact, a lower
    bound on the norm) correct to within a factor 3."

    .. versionadded:: 0.13.0

    References
    ----------
    .. [1] Nicholas J. Higham and Francoise Tisseur (2000),
           "A Block Algorithm for Matrix 1-Norm Estimation,
           with an Application to 1-Norm Pseudospectra."
           SIAM J. Matrix Anal. Appl. Vol. 21, No. 4, pp. 1185-1201.

    .. [2] Awad H. Al-Mohy and Nicholas J. Higham (2009),
           "A new scaling and squaring algorithm for the matrix exponential."
           SIAM J. Matrix Anal. Appl. Vol. 31, No. 3, pp. 970-989.

    
"""
function onenormest(A, p::Int=1, t::Int=2, itmax::Int=5)
    if size(A,1) != size(A,2)
        error("expect square matrix.")
    end

    est, nmults, nresamples = _onenormest_impl(A, adjoint(A), p, t, itmax)
    
    #println("[one norm est] nmults: ", nmults*p)
    return est
end



## _mulp! 
## Y = A^p * X
## p must be positive >=1
function _mulp!(Y::AbstractVecOrMat, A, X::AbstractVecOrMat, p::Int)
    
    # allocate tempo space
    w = similar(X)

    mul!(w,A,X)
    copyto!(Y,w)
    for i in 2:p
        mul!(w,A,Y)
        copyto!(Y,w)
    end
    Y
end

"""
    From Higham and Tisseur:
    "Everything in this section remains valid for complex matrices
    provided that sign(A) is redefined as the matrix (aij / |aij|)
    (and sign(0) = 1) transposes are replaced by conjugate transposes."
"""
function _sign_roundup(X::T) where {T <: Number} 
    if X == 0
        return T(1)
    else
        return X / abs(X)
    end
end



## note this does not check the size of A and B, they need to be the same
@inline function _check_vecs_parallel(a::AbstractVector{T}, b::AbstractVector{T}) where {T <: Number} 
    return abs(dot(a,b)) .≈ length(a)
end


@inline function _check_vec_parallel_in_mat(a::AbstractVector{T}, B::AbstractMatrix{T}) where {T <: Number} 
    return any(abs.(adjoint(a) * B) .≈ length(a))
end

## note this does not check the size of A and B. They need to be the same
## size(A)  == (n,t)
## size(B)  == (n,t)
function _check_cols_parallel(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T <: Number} 
    n = size(A,1)
    for j in 1:size(A,2)
        if !_check_vec_parallel_in_mat(A[:,j],B)
            return false
        end
    end
    return true
end


function _check_need_resample(i, X::AbstractMatrix{T}, Y::AbstractMatrix{T}) where {T <: Number}
    # i-th col of X need to resample if 
    #    1) its parallel for a previous col of X
    # OR 2) its parallel to a col of Y
    n, t = size(X)
    v = X[:,i]
    for j in 1:i-1
        if _check_vecs_parallel(v,X[:,j])
            return true
        end
    end

    for j in 1:t
        if _check_vecs_parallel(v,Y[:,j])
            return true
        end
    end

    return false

end 

@inline function _resample_col(i::Int, X::AbstractMatrix{T}) where {T <: Number}
    X[:,i] = rand((1.,-1.),size(X,1))
end

"""
    Compute a lower bound of the 1-norm of a sparse matrix, [implementation]

    Parameters
    ----------
    A : AbstractMatrix
        A linear operator that can produce matrix products, must be square matrix with (nxn)
    AT : AbstractMatrix
        The transpose of A.
    p : int, optional (default =1)
        The power of A
    t : int, optional
        A positive parameter controlling the tradeoff between
        accuracy versus time and memory usage.
    itmax : int, optional
        Use at most this many iterations.

    Returns
    -------
    est : float
        An underestimate of the 1-norm of the sparse matrix.
 
    nmults : int, optional
        The number of matrix products that were computed.
    nresamples : int, optional
        The number of times a parallel column was observed,
        necessitating a re-randomization of the column.

    [Note] This is algorithm 2.4.

"""
function _onenormest_impl(A, AT, p::Int=1, t::Int=2, itmax::Int=5) 
    if itmax < 2
        error("itmax must be at least 2")
    end
    if t < 1
        error("must be at least one column, t>=1")
    end

    n = size(A,1)

    if t >= n
        error("t must be less than the order of A")
    end

    T = eltype(A)

    nmults::Int = 0
    nresamples::Int = 0

    X = ones(T,n,1)
    if t > 1
        #X = hcat(X,[1.,1.,1.,-1.,1.,-1.,-1.,-1.,1.,-1.]) # testing
        X = hcat(X,rand((1.,-1.),n,t-1))
        ## checking with previous col, to see if resample is needed
        for i in 2:t
            while _check_vecs_parallel(X[:,i],X[:,i-1])
                #need resample i-th column
                X[:,i] = rand((1.,-1.),n)
                nresamples += 1
            end
        end
    end

    # normalize for each column to be unit 1-norm
    X ./= n

    est =0
    est_old = 0
    ind_hist = Vector{Int}(undef,0)
    k::Int = 1
    S = zeros(T,n,t)
    
    ind = nothing
    ind_best = nothing
    
    while true
        #@show X
        Y = similar(X)
        _mulp!(Y,A,X,p) # mulp: Y = A^p * X
        #@show Y

        nmults += 1
        mags = collect( norm(Y[:,j],1) for j in 1:t )
        est, best_j = findmax(mags)
        #println(mags)
        #println(est)
        #println(best_j)
        if (est > est_old) || (k==2)
            if k >= 2
                ind_best = ind[best_j]
            end 
            w = Y[:,best_j]
        end

        # (1)
        if (k >= 2) && (est <= est_old)
            est = est_old
            break
        end
        est_old = est
        S_old = S
        if k > itmax
            break
        end

        S = _sign_roundup.(Y)

        #@show S

        #(2)
        if _check_cols_parallel(S, S_old)
            break
        end

        if t > 1
            for i in 1:t
                while _check_need_resample(i,S,S_old)
                    S[:,i] = rand((1.,-1.),n)
                    nresamples += 1
                end 
            end
        end
        S_old = nothing

        #(3), reuse Y: (note Y is the Z)
        _mulp!(Y,AT,S,p) # mulp: Y = A^p * X

        #display(Y)


        nmults += 1
        Y = abs.(Y)

        h = collect( maximum(Y[j,:]) for j in 1:n )
        #println("h")
        #display(h)


        Y = nothing

        

        #(4) 
        if (k >= 2)  
            if(maximum(h) == h[ind_best])       
                break
            end
        end

        ## sort h in descending order, and re-order ind correspondingly.
        ind = reverse(sortperm(h))[1:t+length(ind_hist)]
        
        h = nothing

        if t > 1
            # (5)
            
            seen = collect(elem in ind_hist for elem in ind)
            # break if the most promising t vec have been visited  
            if all(seen[1:t])
                break
            end
            
            ind = vcat(ind[map(!,seen)], ind[seen])
        end

        fill!(X,0)
        for j in 1:t 
            X[ind[j],j] = 1.
        end

        seen = collect(elem in ind_hist for elem in ind[1:t])
        new_ind = (ind[1:t])[map(!,seen)]
        ind_hist = vcat(ind_hist, new_ind)
        k += 1

    end
    #println(ind_best)
    return est, nmults, nresamples 

end