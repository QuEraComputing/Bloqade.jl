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
Base.:*(n, m::T) where {T <: ThreadedMatrix} = n * m.matrix
Base.pointer(m::T) where {T <: Diagonal} = pointer(m.diag)

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

Adapt.@adapt_structure Hamiltonian

"""
    struct StepHamiltonian

A low-level linear-map object that encodes time-dependent
hamiltonian at time step `t`. This object supports the
linear map interface `mul!(Y, H, X)`.
"""
struct StepHamiltonian{T,FS,TS}
    t::T # clock
    h::Hamiltonian{FS,TS}
end

Base.size(h::StepHamiltonian, idx::Int) = size(h.h, idx)
Base.size(h::StepHamiltonian) = size(h.h)

function to_matrix(h::StepHamiltonian)
    return sum(zip(h.h.fs, h.h.ts)) do (f, t)
        return f(h.t) * t
    end
end

function LinearAlgebra.opnorm(h::StepHamiltonian, p = 2)
    return opnorm(to_matrix(h), p)
end

(h::Hamiltonian)(t::Real) = StepHamiltonian(t, h)


storage_size(x) = sizeof(x)
storage_size(H::T) where {T <: ThreadedMatrix} = storage_size(H.matrix)
function storage_size(h::Hamiltonian)
    return sum(storage_size, h.ts)
end
function storage_size(H::SparseMatrixCSC)
    return sizeof(H.colptr) + sizeof(H.rowval) + sizeof(H.nzval)
end