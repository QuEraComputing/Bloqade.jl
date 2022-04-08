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

# https://github.com/JuliaGPU/CUDA.jl/pull/1470
for SparseMatrixType in [:CuSparseMatrixCSC, :CuSparseMatrixCSR]
    @eval begin
        CUDA.CUSPARSE.$SparseMatrixType(S::Diagonal) = $SparseMatrixType(cu(S))
        CUDA.CUSPARSE.$SparseMatrixType(S::Diagonal{T, <:CuArray}) where T = $SparseMatrixType{T}(S)
        CUDA.CUSPARSE.$SparseMatrixType{Tv}(S::Diagonal{T, <:CuArray}) where {Tv, T} = $SparseMatrixType{Tv, Cint}(S)
        function CUDA.CUSPARSE.$SparseMatrixType{Tv, Ti}(S::Diagonal{T, <:CuArray}) where {Tv, Ti, T}
            m = size(S, 1)
            return $SparseMatrixType{Tv, Ti}(CuVector(1:(m+1)), CuVector(1:m), Tv.(S.diag), (m, m))
        end
    end
end
