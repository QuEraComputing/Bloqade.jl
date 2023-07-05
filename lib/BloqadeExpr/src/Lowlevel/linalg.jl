#**************************************************************
#  Here, one can find the linalg API for  
#  1. Hamiltonian/SumOfLinop
#  2. Binding of standard LinearAlgebra functions 
#     to the backend of choice for ThreadedMatrix
#**************************************************************


## mul!
#--------------------------------
function LinearAlgebra.mul!(C::AbstractVecOrMat, A::SumOfLinop, B::AbstractVecOrMat)
    fill!(C, zero(eltype(C)))
    #print(eltype(C))
    for (f, term) in zip(A.fvals, A.ts)
        #println("terms:", typeof(term), " ", typeof(f))
        mul!(C, term, B, f, one(f))
    end
    return C
end

## additionals, maybe we don't need this.
function Base.:*(a::Number,b::SumOfLinop)
    return SumOfLinop{RegularLinop}(b.fvals .* a, b.ts)
end

function Base.:*(a::Complex,b::SumOfLinop{T}) where {T}
    if real(a) ≈ 0
        return SumOfLinop{anti_type(T)}(b.fvals .* a, b.ts)
    elseif imag(a) ≈ 0
        return SumOfLinop{T}(b.fvals .* a, b.ts)
    else
        return SumOfLinop{RegularLinop}(b.fvals .* a, b.ts)
    end
end

function Base.:*(a::Real,b::SumOfLinop{T}) where {T}
    return SumOfLinop{T}(b.fvals .* a, b.ts)
end

Base.:*(n, m::T) where {T <: ThreadedMatrix} = n * m.matrix


#=
function Base.:+(a::SumOfLinop, b::SumOfLinop)
    if !(a.ts === b.ts)
        error("two SumOfLinop must share the same static terms ")
    end
    return SumOfLinop(a.fvals + b.fvals, a.ts)
end

function Base.:-(a::SumOfLinop, b::SumOfLinop)
    if !(a.ts === b.ts)
        error("two SumOfLinop must share the same static terms ")
    end
    return SumOfLinop(a.fvals - b.fvals, a.ts)
end
=#


# if BloqadeExpr was the backend of choice, then ThreadedMatrix will have the SparseMatrixCSC type
# in which case we just dispatch to standard mul!
LinearAlgebra.mul!(C, A::ThreadedMatrix{<:SparseMatrixCSC}, B, α, β) = SparseArrays.mul!(C, A.matrix, B, α, β) 

# if PMCSR/TSCSR was the backend of choice, then direct to bmul!'s below
LinearAlgebra.mul!(C, A::ThreadedMatrix, B, α, β) = bmul!(C, A.matrix, B, α, β)

##-------------------------------- mul!





## opnorm()
# --------------------------------
function LinearAlgebra.opnorm(h::SumOfLinop, p = 2)
    return opnorm(to_matrix(h), p)
end

##---------------------------------


## tr()
# --------------------------------
function LinearAlgebra.tr(A::SumOfLinop)
    return sum(zip(A.fvals, A.ts)) do (f, t)
        return f * tr(t)
    end
end

# [TODO] parallel trace (btrace)
# if BloqadeExpr was the backend of choice, then ThreadedMatrix will have the SparseMatrixCSC type
# in which case we just dispatch to single threaded tr
##LinearAlgebra.tr(C, A::ThreadedMatrix{<:SparseMatrixCSC}, B, α, β) = trace(A.matrix) 

# if PMCSR/TSCSR was the backend of choice, then direct to bmul!'s below
##LinearAlgebra.tr(C, A::ThreadedMatrix, B, α, β) = btrace(A.matrix)

## [NOTE] currently only single thread for trace is implemented, so we use the exposed LinearAlgebra.tr() 
#  that bind to trace()
LinearAlgebra.tr(A::ThreadedMatrix) = tr(A.matrix)

##--------------------------------  tr()

## check if is hermitian. 
# --------------------------------
LinearAlgebra.ishermitian(A::SumOfLinop{<: LinearAlgebra.Hermitian}) = true
LinearAlgebra.ishermitian(A::SumOfLinop) = false

isskewhermitian(A::SumOfLinop{<: SkewHermitian}) = true
isskewhermitian(A::SumOfLinop) = false


## adjoint() 
function LinearAlgebra.adjoint(A::SumOfLinop{<:LinearAlgebra.Hermitian})
    return A
end
function LinearAlgebra.adjoint(A::SumOfLinop{<:SkewHermitian})
    return SumOfLinop{SkewHermitian}(A.fvals.*(-1), A.ts)
end
function LinearAlgebra.adjoint(A::SumOfLinop{OPTYPE}) where {OPTYPE}
    return SumOfLinop{OPTYPE}(conj.(A.fvals), map(adjoint,A.ts))
end

## add constant identity term into SumOfLinop
# [NOTE] this does not check the type consistency of c w.r.t. A.fvals. 
function add_I(A,c::Number)
    Iop = LinearAlgebra.I(size(A,1))
    return A + c*Iop
end

function add_I(A::SumOfLinop, c::Number)
    Iop = LinearAlgebra.I(size(A,1))

    if nthreads() > 1
        return SumOfLinop{RegularLinop}((A.fvals...,c),(A.ts...,ThreadedMatrix(Iop)))
    else
        return SumOfLinop{RegularLinop}((A.fvals...,c),(A.ts...,Iop))
    end

end

function add_I(A::SumOfLinop{<:LinearAlgebra.Hermitian}, c::Real)
    # check backend:
    Iop = LinearAlgebra.I(size(A,1))

    if nthreads() > 1
        return SumOfLinop{LinearAlgebra.Hermitian}((A.fvals...,c),(A.ts...,ThreadedMatrix(Iop)))
    else
        return SumOfLinop{LinearAlgebra.Hermitian}((A.fvals...,c),(A.ts...,Iop))
    end

end
function add_I(A::SumOfLinop{<:LinearAlgebra.Hermitian}, c::Complex)
    # check backend:
    Iop = LinearAlgebra.I(size(A,1))

    OPTYPE=RegularLinop
    if imag(c) ≈ 0
        OPTYPE = LinearAlgebra.Hermitian
    end

    if nthreads() > 1
        return SumOfLinop{OPTYPE}((A.fvals...,c),(A.ts...,ThreadedMatrix(Iop)))
    else
        return SumOfLinop{OPTYPE}((A.fvals...,c),(A.ts...,Iop))
    end

end

function add_I(A::SumOfLinop{<:SkewHermitian}, c::Real)
    # check backend:
    Iop = LinearAlgebra.I(size(A,1))

    if nthreads() > 1
        return SumOfLinop{RegularLinop}((A.fvals...,c),(A.ts...,ThreadedMatrix(Iop)))
    else
        return SumOfLinop{RegularLinop}((A.fvals...,c),(A.ts...,Iop))
    end

end
function add_I(A::SumOfLinop{<:SkewHermitian}, c::Complex)
    # check backend:
    Iop = LinearAlgebra.I(size(A,1))

    OPTYPE=RegularLinop
    if real(c) ≈ 0
        OPTYPE = SkewHermitian
    end

    if nthreads() > 1
        return SumOfLinop{OPTYPE}((A.fvals...,c),(A.ts...,ThreadedMatrix(Iop)))
    else
        return SumOfLinop{OPTYPE}((A.fvals...,c),(A.ts...,Iop))
    end

end