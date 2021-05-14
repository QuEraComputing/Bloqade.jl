struct ExpMCache{T}
    A2::Matrix{T}
    P::Matrix{T}
    U::Matrix{T}
    V::Matrix{T}
    T::Matrix{T}
end

function ExpMCache(A::StridedMatrix{T}) where T
    n = LinearAlgebra.checksquare(A)
    A2   = Matrix{T}(undef, n, n)
    P    = Matrix{T}(undef, n, n)
    U    = Matrix{T}(undef, n, n)
    V    = Matrix{T}(undef, n, n)
    temp = Matrix{T}(undef, n, n)
    return ExpMCache(A2, P, U, V, temp)
end

## Destructive matrix exponential using algorithm from Higham, 2008,
## "Functions of Matrices: Theory and Computation", SIAM
##
## Non-allocating version of `LinearAlgebra.exp!`. Modifies `A` to
## become (approximately) `exp(A)`.
function exp!(A::StridedMatrix{T}, cache::ExpMCache=ExpMCache(A)) where T <: LinearAlgebra.BlasFloat
    X = A
    n = LinearAlgebra.checksquare(A)
    # if ishermitian(A)
        # return copytri!(parent(exp(Hermitian(A))), 'U', true)
    # end

    A2, P, U, V, temp = cache.A2, cache.P, cache.U, cache.V, cache.T
    fill!(P, zero(T));
    for k in 1:n
        P[k, k] = one(T)
    end

    ilo, ihi, scale = LAPACK.gebal!('B', A)    # modifies A
    nA = opnorm(A, 1)
    ## For sufficiently small nA, use lower order Padé-Approximations
    if (nA <= 2.1)
        if nA > 0.95
            C = T[17643225600.,8821612800.,2075673600.,302702400.,
                     30270240.,   2162160.,    110880.,     3960.,
                           90.,         1.]
        elseif nA > 0.25
            C = T[17297280.,8648640.,1995840.,277200.,
                     25200.,   1512.,     56.,     1.]
        elseif nA > 0.015
            C = T[30240.,15120.,3360.,
                    420.,   30.,   1.]
        else
            C = T[120.,60.,12.,1.]
        end
        mul!(A2, A, A)
        @. U = C[2] * P
        @. V = C[1] * P
        for k in 1:(div(size(C, 1), 2) - 1)
            k2 = 2 * k
            mul!(temp, P, A2); P, temp = temp, P # equivalent to P *= A2
            @. U += C[k2 + 2] * P
            @. V += C[k2 + 1] * P
        end
        mul!(temp, A, U); U, temp = temp, U # equivalent to U = A * U
        @. X = V + U
        @. temp = V - U
        LAPACK.gesv!(temp, X)
    else
        s  = log2(nA/5.4)               # power of 2 later reversed by squaring
        if s > 0
            si = ceil(Int,s)
            A ./= convert(T,2^si)
        end
        C  = T[64764752532480000.,32382376266240000.,7771770303897600.,
                1187353796428800.,  129060195264000.,  10559470521600.,
                    670442572800.,      33522128640.,      1323241920.,
                        40840800.,           960960.,           16380.,
                             182.,                1.]
        mul!(A2, A, A)
        @. U = C[2] * P
        @. V = C[1] * P
        for k in 1:6
            k2 = 2 * k
            mul!(temp, P, A2); P, temp = temp, P # equivalent to P *= A2
            @. U += C[k2 + 2] * P
            @. V += C[k2 + 1] * P
        end
        mul!(temp, A, U); U, temp = temp, U # equivalent to U = A * U
        @. X = V + U
        @. temp = V - U
        LAPACK.gesv!(temp, X)

        if s > 0            # squaring to reverse dividing by power of 2
            for t=1:si
                mul!(temp, X, X)
                X .= temp
            end
        end
    end

    # Undo the balancing
    for j = ilo:ihi
        scj = scale[j]
        for i = 1:n
            X[j,i] *= scj
        end
        for i = 1:n
            X[i,j] /= scj
        end
    end

    if ilo > 1       # apply lower permutations in reverse order
        for j in (ilo-1):-1:1; LinearAlgebra.rcswap!(j, Int(scale[j]), X) end
    end
    if ihi < n       # apply upper permutations in forward order
        for j in (ihi+1):n;    LinearAlgebra.rcswap!(j, Int(scale[j]), X) end
    end
    return X
end
