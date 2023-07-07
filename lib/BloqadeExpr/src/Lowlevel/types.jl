# a linear map for low-level hamiltonian representation
const backend = @load_preference("backend", "BloqadeExpr")

struct ThreadedMatrix{M <: AbstractMatrix}
    matrix::M

    # use Inner Constructor to avoid infinite loop of calling constructor on itself
    ThreadedMatrix(m::T) where {T<:AbstractMatrix} = new{typeof(m)}(m)
    function ThreadedMatrix(m::SparseMatrixCSC)
        @static if backend == "ParallelMergeCSR" # should be conjugate transpose
            transformed_matrix =  m |> conj! |> transpose
        elseif backend == "ThreadedSparseCSR" # should be conjugate transpose, then turned into SparseMatrixCSR
            transformed_matrix = m |> conj! |> transpose |> SparseMatrixCSR 
        elseif backend == "BloqadeExpr"
            transformed_matrix = m
        else
            throw(ArgumentError("The backend selected is not supported."))
        end

        return new{typeof(transformed_matrix)}(transformed_matrix)
    end

end

Base.size(m::ThreadedMatrix) = size(m.matrix)
Base.size(m::ThreadedMatrix, i) = size(m.matrix)[i]
Base.pointer(m::T) where {T <: Diagonal} = pointer(m.diag)


precision_type(m::T) where {T <: Number} = real(typeof(m))
precision_type(m::T) where {T <: Diagonal} = real(eltype(m))
precision_type(m::T) where {T <: PermMatrix} = real(eltype(m))
precision_type(m::T) where {T <: SparseMatrixCSR} = real(eltype(m))
precision_type(m::T) where {T <: SparseMatrixCSC} = real(eltype(m))
precision_type(m::T) where {T <: ThreadedMatrix} = real(eltype(m.matrix))


"""
    struct Hamiltonian

`Hamiltonian` stores the dynamic prefactors of each term.
The actual hamiltonian is the sum of `f_i(t) * t_i` where
`f_i` and `t_i` are entries of `fs` and `ts`.
"""
struct Hamiltonian{FS<:Tuple,TS<:Tuple}
    fs::FS # prefactor of each term
    ts::TS # const linear map of each term

    function Hamiltonian(fs, ts)
        all(x -> size(x) == size(ts[1]), ts) || throw(ArgumentError("matrix term should have the same size"))
        return new{typeof(fs),typeof(ts)}(fs, ts)
    end
end

Base.size(h::Hamiltonian) = size(h.ts[1])
Base.size(h::Hamiltonian, idx::Int) = size(h.ts[1], idx)
function precision_type(h::Hamiltonian) 
    tp = unique(precision_type.(h.ts))
    return Union{tp...}
end
function highest_type(h::Hamiltonian)
    tp = unique(eltype.(h.ts))
    return promote_type(tp...)
end


Base.eltype(h::Hamiltonian) = highest_type(h)


Adapt.@adapt_structure Hamiltonian




abstract type RegularLinop end
abstract type SkewHermitian end

anti_type(::Type{LinearAlgebra.Hermitian}) = SkewHermitian
anti_type(::Type{SkewHermitian}) = LinearAlgebra.Hermitian
anti_type(::Type{RegularLinop}) = RegularLinop


"""
    struct SumOfLinop
A low-level linear-map object that explicitly evaluate time dependent 
coefficients at given time `t` fvals = fs(t) of Hamiltonian. 

This object supports the linear map interface `mul!(Y, H, X)`.
"""

struct SumOfLinop{OPTYPE, VS,TS}
    fvals::VS
    ts::TS
    function SumOfLinop{OPTYPE}(fvals::VS, ts::TS) where {OPTYPE, VS, TS}
        return new{OPTYPE,VS,TS}(fvals, ts)
    end
end

Base.size(h::SumOfLinop, idx::Int) = size(h.ts[1], idx)
Base.size(h::SumOfLinop) = size(h.ts[1])
function precision_type(h::SumOfLinop)
    tp = unique(precision_type.(h.ts))
    tp2 = unique(precision_type.(h.fvals))
    tp = unique((tp...,tp2...))
    return Union{tp...}
end
function highest_type(h::SumOfLinop)
    tp = unique(eltype.(h.ts))
    tp2 = unique(typeof.(h.fvals))
    return promote_type(tp...,tp2...)
end
Base.eltype(h::SumOfLinop) = highest_type(h)

function to_matrix(h::SumOfLinop)
    return sum(zip(h.fvals, h.ts)) do (f, t)
        return f * t
    end
end

function _getf(h::Hamiltonian,t) 
    return collect(map(h.fs) do f 
        return f(t)
    end
    )
end


## lowering by Hamiltonian, so its Hermitian type
(h::Hamiltonian)(t::Real) = SumOfLinop{LinearAlgebra.Hermitian}(_getf(h,t), h.ts)




storage_size(x) = sizeof(x)
storage_size(H::T) where {T <: ThreadedMatrix} = storage_size(H.matrix)
function storage_size(h::Hamiltonian)
    return sum(storage_size, h.ts)
end
function storage_size(H::SparseMatrixCSC)
    return sizeof(H.colptr) + sizeof(H.rowval) + sizeof(H.nzval)
end
