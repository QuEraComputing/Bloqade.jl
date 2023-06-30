using Revise
using BloqadeKrylov
using Bloqade
using LinearAlgebra
using BloqadeLattices 


function writeSparseCSC(A,fname)

    open("$fname.size","w") do f
        write(f,A.m)
        write(f,A.n)
    
    end
    
    open("$fname.colptr","w") do f
        write(f,A.colptr)
    end
    open("$fname.rowval","w") do f
        write(f,A.rowval)
    end
    open("$fname.nzval","w") do f
        write(f,A.nzval)
    end

end

## wrapping Hamiltonian with counter that count the number of mul! calls
mutable struct WHam
    h
    counter::Vector{Int}
    function WHam(h, counter= Vector{Int}([0]))
        return new(h,counter)
    end
end

Base.size(M::WHam) = size(M.h)
Base.size(M::WHam,i::Int) = size(M.h,i)
#opnorm(M::WHam) = opnorm(M.h)
LinearAlgebra.opnorm(M::WHam,p=2) = opnorm(M.h,p)
Base.eltype(M::WHam) = eltype(M.h)


function count_mul!(a::AbstractVecOrMat, M::WHam, x::AbstractVecOrMat)
    M.counter[1] += 1
    mul!(a, M.h, x)
end
LinearAlgebra.mul!(a::AbstractVecOrMat,M::WHam,x::AbstractVecOrMat) = count_mul!(a,M,x)
LinearAlgebra.tr(M::WHam) = tr(M.h)
BloqadeExpr.add_I(M::WHam, c::Number) =  WHam(BloqadeExpr.add_I(M.h, c), M.counter)
Base.:*(c::Number, M::WHam) = WHam(c*M.h, M.counter)
LinearAlgebra.adjoint(M::WHam) = WHam(adjoint(M.h), M.counter)
##===================================


# reg inplace modify
function testing_expmv!(H::BloqadeExpr.Hamiltonian, tspan, reg)
    n_mul = Vector{Int}(undef,0)
    for t in tspan
        Ht = WHam(H(t))
    
        BloqadeKrylov.expmv!((dt)*im, Ht, statevec(reg), tol=1.0e-7)
    
        #println(norm(statevec(reg)))
        #println(Ht.counter[1])
        push!(n_mul, Ht.counter[1])
    end
    return n_mul
end

function testing_expm_multiply!(H::BloqadeExpr.Hamiltonian, tspan, reg)
    n_mul_impl = Vector{Int}(undef,0)
    n_mul_sm = Vector{Int}(undef,0)
    print(tspan)
    for t in tspan
        Ht = H(t)
        Ht = WHam(Ht)
    

        #BloqadeKrylov.expm_multiply!((dt)*im, Ht, statevec(reg1))
        #^> unwrap it to get counting for mul from get optimal sm overhead
        ti = (dt)


        m_star, s, μ, A_1norm, As = BloqadeKrylov.get_optimal_sm(ti,im*Ht)
        push!(n_mul_sm, Ht.counter[1])
        Ht.counter[1] = 0

        v = statevec(reg)
        w = similar(v)


        #assuming its ComplexF64
        BloqadeKrylov._expm_multiply_impl!(w, ti, As, v, μ, s, m_star,1.0e-7)



        copyto!(v,w)
        
        
        #println(norm(statevec(reg)))
        #println(Ht.counter[1])
        push!(n_mul_impl, Ht.counter[1])
    end
    return n_mul_sm, n_mul_impl
end
