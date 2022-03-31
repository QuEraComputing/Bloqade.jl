# NOTE: copied from Expokit
# TODO: fully switch to ExponentialUtilities after fixing complex number
function default_anorm(A)
    try
        opnorm(A, Inf)
    catch err
        if err isa MethodError
            @warn "opnorm($(typeof(A)), Inf) is not defined, fall back to using `anorm = 1.0`.
To suppress this warning, please specify the `anorm` keyword manually."
            1.0
        else
            throw(err)
        end
    end
end

"""
    expmv(t, A, vec; [tol], [m], [norm], [anorm])
Calculate matrix exponential acting on some vector, ``w = e^{tA}v``,
using the Krylov subspace approximation.
See R.B. Sidje, ACM Trans. Math. Softw., 24(1):130-156, 1998
and http://www.maths.uq.edu.au/expokit
"""
function expmv( t::Number,
                   A, vec::Vector{T};
                   tol::Real=1e-7,
                   m::Int=min(30, size(A, 1)),
                   norm=LinearAlgebra.norm, anorm=default_anorm(A)) where {T}
    result = convert(Vector{promote_type(eltype(A), T, typeof(t))}, copy(vec))
    expmv!(t, A, result; tol=tol, m=m, norm=norm, anorm=anorm)
    return result
end

expmv!( t::Number,
    A, vec::Vector{T};
    tol::Real=1e-7,
    m::Int=min(30,size(A,1)),
    norm=LinearAlgebra.norm, anorm=default_anorm(A)) where {T} = expmv!(vec, t, A, vec; tol=tol, m=m, norm=norm, anorm=anorm)

function expmv!(w::Vector{T}, t::Number, A, vec::Vector{T};
                   tol::Real=1e-7, m::Int=min(30,size(A,1)), norm=LinearAlgebra.norm, anorm=default_anorm(A)) where {T}
    if size(vec,1) != size(A,2)
        error("dimension mismatch")
    end
    # safety factors
    gamma = 0.9
    delta = 1.2
    btol = 1e-7     # tolerance for "happy-breakdown"
    maxiter = 10    # max number of time-step refinements
    rndoff= anorm*eps()
    # estimate first time-step and round to two significant digits
    beta = norm(vec)
    r = 1/m
    fact = (((m+1)/exp(1.0))^(m+1))*sqrt(2.0*pi*(m+1))
    tau = (1.0/anorm)*((fact*tol)/(4.0*beta*anorm))^r
    tau = round(tau, sigdigits=2)
    # storage for Krylov subspace vectors
    vm = Array{typeof(w)}(undef,m+1)
    for i=1:m+1
        vm[i]=similar(w)
    end
    hm = zeros(T,m+2,m+2)
    tf = abs(t)
    tsgn = sign(t)
    tk = zero(tf)
    copyto!(w, vec)
    p = similar(w)
    mx = m
    while tk < tf
        tau = min(tf-tk, tau)
        # Arnoldi procedure
        # vm[1] = v/beta
        rmul!(copyto!(vm[1],w),1/beta)
        mx = m
        for j=1:m
            # p[:] = A*vm[j]
            mul!(p, A, vm[j])
            for i=1:j
                hm[i,j] = dot(vm[i], p)
                # p[:] = p - hm[i,j]*vm[i]
                p = axpy!(-hm[i,j], vm[i], p)
            end
            s = norm(p)
            if s < btol # happy-breakdown
                tau = tf - tk
                err_loc = btol
                # F = expm(tsgn*tau*hm[1:j,1:j])
                # F = expm!(scale(tsgn*tau,view(hm,1:j,1:j)))
                L = tsgn*tau*view(hm,1:j,1:j)
                F = ExponentialUtilities.exponential!(L)
                fill!(w, zero(T))
                for k=1:j
                    # w[:] = w + beta*vm[k]*F[k,1]
                    w = axpy!(beta*F[k,1], vm[k], w)
                end
                mx = j
                break
            end
            hm[j+1,j] = s
            # vm[j+1] = p/hm[j+1,j]
            rmul!(copyto!(vm[j+1],p),1/hm[j+1,j])
        end
        hm[m+2,m+1] = one(T)
        (mx != m) || (avnorm = norm(mul!(p,A,vm[m+1])))
        # propagate using adaptive step size
        iter = 1
        while (iter < maxiter) && (mx == m)
            # F = expm(tsgn*tau*hm)
            # F = expm!(scale(tsgn*tau,hm))
            F = ExponentialUtilities.exponential!(tsgn*tau*hm)
            # local error estimation
            err1 = abs( beta*F[m+1,1] )
            err2 = abs( beta*F[m+2,1] * avnorm )
            if err1 > 10*err2	# err1 >> err2
                err_loc = err2
                r = 1/m
            elseif err1 > err2
                err_loc = (err1*err2)/(err1-err2)
                r = 1/m
            else
                err_loc = err1
                r = 1/(m-1)
            end
            # time-step sufficient?
            if err_loc <= delta * tau * (tau*tol/err_loc)^r
                fill!(w, zero(T))
                for k=1:m+1
                    # w[:] = w + beta*vm[k]*F[k,1]
                    w = axpy!(beta*F[k,1], vm[k], w)
                end
                break
            end
            tau = gamma * tau * (tau*tol/err_loc)^r # estimate new time-step
            tau = round(tau, sigdigits=2) # round to 2 signiﬁcant digits
                                 # to prevent numerical noise
            iter = iter + 1
        end
        if iter == maxiter
            # TODO, here an exception should be thrown, but which?
            error("Number of iterations exceeded $(maxiter). Requested tolerance might be to high.")
        end
        beta = norm(w)
        tk = tk + tau
        tau = gamma * tau * (tau*tol/err_loc)^r # estimate new time-step
        tau = round(tau, sigdigits=2) # round to 2 signiﬁcant digits
                             # to prevent numerical noise
        err_loc = max(err_loc,rndoff)
        fill!(hm, zero(T))
    end
    return w
end # expmv!
