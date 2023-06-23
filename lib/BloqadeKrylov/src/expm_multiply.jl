function expm_multiply(t::Number
                       A
                       v::AbstractVector{T},
                       tol::Real = nothing
                       )
    if size(v, 1) != size(A, 2)
        error("dimension mismatch")
    end


    if tol === nothing
        # infer tol base on the size of T
        size_T = sizeof(T)
        if size_T == 8
            tol = 2^(-53)
        elseif size_T == 4
            tol = 2^(-24)
        else
            error("Unsupported type of vec{T}")
        end
    end
    n = size(A,1)
    n0 = 1

    # 1) shift the A:
    traceA = LinearAlgebra.tr(A)
    mu = traceA/n

    As = A - mu*LinearAlgebra.I(n)

    # 2) get the optimal s and m_star
    
    # 2) here lies the impl call:
    F = deepcopy(v)
    _expm_multiply_impl!(F, t, As, v, mu, s, m_star, tol)
end 
                        



"""
    _expm_multiply_impl(t, A, vec; [traceA])
    Calculate matrix exponential acting on some vector, ``w = e^{tA}v``,
    using the Krylov subspace approximation.
"""
function _expm_multiply_impl!(F::AbstractVector{T}
                             t::Number,
                             A, 
                             v::AbstractVector{T},
                             mu::Real, #shift parameter
                             s::Int,
                             m_star::Int,
                             tol::Real,
                            ) where {T}


   

    mu = traceA/n
     
    Ashift = A - mu*LinearAlgebra.I(n)



end