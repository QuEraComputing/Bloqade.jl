@static if !hasmethod(+, Tuple{Diagonal{T,<:CuArray}, CuSparseMatrixCSC{T}} where T)
# https://github.com/JuliaGPU/CUDA.jl/issues/1469

for SparseMatrixType in [:CuSparseMatrixCSC, :CuSparseMatrixCSR]
    @eval begin
        function Base.:+(lhs::Diagonal{T,<:CuArray}, rhs::$SparseMatrixType{T}) where {T}
            return $SparseMatrixType(lhs) + rhs
        end
        function Base.:+(lhs::$SparseMatrixType{T}, rhs::Diagonal{T,<:CuArray}) where {T}
            return lhs + $SparseMatrixType(rhs)
        end
    end
end

end

function CuYao.cu(reg::SubspaceArrayReg{D}) where {D}
    println("here")
    natoms = reg.natoms
    new_state = CuArray(reg.state)
    new_subspace = Subspace(
        reg.subspace.nqubits,
        reg.subspace.map,
        CuArray(reg.subspace.subspace_v)
    )
    return SubspaceArrayReg{D}(new_state,new_subspace)
end