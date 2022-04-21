# https://github.com/JuliaGPU/CUDA.jl/issues/1469

for SparseMatrixType in [:CuSparseMatrixCSC, :CuSparseMatrixCSR]
    @eval begin
        function Base.:+(lhs::Diagonal{T, <:CuArray}, rhs::$SparseMatrixType{T}) where T
            return $SparseMatrixType(lhs) + rhs
        end
        function Base.:+(lhs::$SparseMatrixType{T}, rhs::Diagonal{T, <:CuArray}) where T
            return lhs + $SparseMatrixType(rhs)
        end 
    end
end
