"""
    struct ValHamiltonian

A low-level linear-map object that explicitly evaluate time dependent 
coefficients at given time `t` fvals = fs(t) of StepHamiltonian. 

This object supports the linear map interface `mul!(Y, H, X)`.

"""
struct ValHamiltonian{VS,FS,TS}
    fvals::VS
    h::Hamiltonian{FS,TS}
end

Base.size(h::ValHamiltonian, idx::Int) = size(h.h, idx)
Base.size(h::ValHamiltonian) = size(h.h)

function to_matrix(h::ValHamiltonian)
    return sum(zip(h.fvals, h.h.ts)) do (f, t)
        return f * t
    end
end

function LinearAlgebra.opnorm(h::ValHamiltonian, p = 2)
    return opnorm(to_matrix(h), p)
end

function get_f(h::StepHamiltonian) 
    return collect(map(h.h.fs) do f 
        return f(h.t)
    end
    )
end

## convert StepHamiltonian to ValHamiltonian 
function Val(h::StepHamiltonian)
    # get values of fs:
    return ValHamiltonian(get_f(h),h.h)
end


function LinearAlgebra.mul!(C::AbstractVecOrMat, A::ValHamiltonian, B::AbstractVecOrMat)
    fill!(C, zero(eltype(C)))
    for (f, term) in zip(A.fvals, A.h.ts)
        mul!(C, term, B, f, one(f))
    end
    return C
end

function Base.:*(a::Number, b::ValHamiltonian)
    return ValHamiltonian(b.fvals .* a, b.h)
end

function Base.:+(a::ValHamiltonian, b::ValHamiltonian)
    if !(a === b)
        error("two ValHamiltonian must share the same static terms ")
    end
    return ValHamiltonian(a.fvals + b.fvals, a.h)
end

function Base.:-(a::ValHamiltonian, b::ValHamiltonian)
    if !(a === b)
        error("two ValHamiltonian must share the same static terms ")
    end
    return ValHamiltonian(a.fvals - b.fvals, a.h)
end

