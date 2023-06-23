
function onenormest_explicit(A::AbstractMatrix, p::Int=1)
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


#=
function onenormest(A::AbstractMatrix, t::Int=2, itmax::Int=5)
    n = size(A,2)
    if t >= n
        error("t must be less than n, t >= n is not implmemented yet.")
    else 
        _onenormest_impl(A, adjoint(A), t, itmax)
    end
end

function _onenormest_impl(A::AbstractMatrix, AT::AbstractMatrix, t::Int, itmax::Int)
    return 0
end
=#