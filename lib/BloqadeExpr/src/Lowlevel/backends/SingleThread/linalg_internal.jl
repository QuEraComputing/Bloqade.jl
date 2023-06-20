
## M::Diagonal already have impl from stdlib LinearAlgebra
# LinearAlgebra.tr(M::Diagonal)  
LinearAlgebra.tr(M::PermMatrix) = trace(M)
LinearAlgebra.tr(M::SparseMatrixCSC) = trace(M)
LinearAlgebra.tr(M::SparseMatrixCSR) = trace(M)

function trace(A::PermMatrix)
    A_vals = A.vals
    A_perm = A.perm
    out = 0
    for (i, perm) in enumerate(A_perm)
        if i == perm
            out += A_vals[i]
        end
    end
    return out
end

function trace(A::SparseMatrixCSC)
    nzvals = A.nzval
    colptr = A.colptr
    rowval = A.rowval

    out = 0
    for j in 1:(length(colptr)-1)
        if colptr[j] == colptr[j+1]
            continue
        end
        for (n,i) in enumerate(rowval[colptr[j]:colptr[j+1]-1])
            if j == i
                out += nzvals[colptr[j]+n-1]     
            end
        end 
    end
    return out
end

function trace(A::SparseMatrixCSR)
    
    nzvals = A.nzval
    colval = A.colval
    rowptr = A.rowptr

    out = 0
    for i in 1:(length(rowptr)-1)
        if rowptr[i] == rowptr[i+1]
            continue
        end
        for (n,j) in enumerate(colval[rowptr[i]:rowptr[i+1]-1])
            if j == i
                out += nzvals[rowptr[i]+n-1]     
            end
        end 
    end
    
    
    return out
end

# dealing with ParallelMergeCSR backend
trace(A::Transpose{<:Any, <:SparseMatrixCSC}) = trace(transpose(A)) 
