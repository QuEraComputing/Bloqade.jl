struct ExpokitOptions{T <: Real, F <: AbstractMatrix, AvT <: AbstractVector, KsT <: KrylovSubspace, E}
    tol::T
    krylov_dim::Int
    max_iter::Int

    # cache
    expH::F
    Av::AvT
    Ks::KsT
    expFcache::E
    expHcache::E
end

function promote_expv(t, A, b)
    promote_type(eltype(A), eltype(b), typeof(t))
end

function ExpokitOptions(t, A, b; tol::Real=1e-7, krylov_dim=min(30, size(A, 1)), max_iter::Int=10)
    Ks = KrylovSubspace(t, A, b, krylov_dim, tol)
    T = promote_expv(t, A, b)
    expH = similar(Ks.H, T)
    Av = similar(Ks.V, (length(b), ))
    return ExpokitOptions(tol, krylov_dim, max_iter, expH, Av, Ks, ExpMCache(expH), ExpMCache(Ks.Hm))
end

struct ExpokitIteration{T, S, WT, MT, VT, Options <: ExpokitOptions}
    tsgn::S
    tf::T
    w::WT
    A::MT
    b::VT
    options::Options
end

function ExpokitIteration(w, t, A, b, options)
    return ExpokitIteration(sign(t), abs(t), w, A, b, options)
end

function Base.iterate(it::ExpokitIteration)
    m, tol = it.options.krylov_dim, it.options.tol
    anorm = opnorm(it.A, Inf); β = norm(it.b)
    tk = zero(eltype(it.tf)); r = 1/m

    fact = (((m+1)/ℯ)^(m+1))*sqrt(2π*(m+1))
    τ0 = (1/anorm)*((fact*tol)/(4β*anorm))^r
    τ0 = round(τ0, sigdigits=2)
    τ0 = min(it.tf - tk, τ0)

    return iterate(it, (τ0, tk, it.b, β))
end

function Base.iterate(it::ExpokitIteration, (τ, tk, b, β))
    w, A, Ks, options = it.w, it.A, it.options.Ks, it.options
    tf, tsgn = it.tf, it.tsgn
    m = it.options.krylov_dim

    tk < tf || return

    krylov_subspace!(Ks, A, b, β)
    if Ks.krylov_dim < options.krylov_dim # happy-breakdown
        happy_breakdown!(w, tk, tf, tsgn, options)
        return
    end

    τ = min(tf-tk, τ)
    @inbounds avnorm = norm(mul!(options.Av, A, Ks.v[end]))

    for _ in 1:options.max_iter
        copyto!(options.expH, Ks.H)
        F = exp!(lmul!(tsgn * τ, options.expH), options.expFcache)
        err, r = local_error_estimation(options.krylov_dim, β, avnorm, F)
        τ_new = estimate_timestep(τ, err, r, options.tol)
        if err ≤ error_upperbound(τ, err, r, options.tol)
            lmul!(β, mul!(w, Ks.V, view(F, 1:m+1, 1)))
            return w, (τ_new, tk + τ, w, norm(w))
        end
        τ = τ_new
    end
    error("maximum iteration exceeded ($max_iter), requested tolerance might be too high.")
end

function expv(t::Number, A, b; kwargs...)
    expv!(similar(b, promote_expv(t, A, b)), t, A, b; kwargs...)
end

function expv!(w, t::Number, A, b; tol::Real=1e-7, krylov_dim=min(30, size(A, 1)), max_iter=10)
    options = ExpokitOptions(t, A, b;tol, krylov_dim, max_iter)
    expv!(w, t, A, b, options)
    return w
end

function expv!(w, t::Number, A, b, options::ExpokitOptions)
    for _ in ExpokitIteration(w, t, A, b, options)
    end
    return w
end

function happy_breakdown!(w, tk::Real, tf::Real, tsgn, options::ExpokitOptions)
    Ks = options.Ks
    tau = tf - tk
    m = Ks.krylov_dim
    # NOTE: the found subspace is smaller than we allocate
    # and since it is exact always use a view of it
    F = exp!(lmul!(tsgn * tau, Ks.Hm), options.expHcache)
    # just use the memory of H, we don't need it anymore
    mul!(w, view(Ks.V, :, 1:30), view(F, :, 1))
    return w
end

function local_error_estimation(m::Int, β::Real, avnorm::Real, F)
    @inbounds begin
        err1 = β * abs(F[m+1, 1])
        err2 = β * abs(F[m+2, 1]) * avnorm
    end

    if err1 > 10err2
        err = err2
        r = 1/m
    elseif err1 > err2
        err = (err1 * err2)/(err1 - err2)
        r = 1/m
    else
        err = err1
        r = 1/(m-1)
    end
    return err, r
end

function estimate_timestep(τ, err, r, tol=1e-7, γ=0.9)
    return round(γ * τ * (τ*tol/err)^r, sigdigits=2)
end

function error_upperbound(τ, err, r, tol=1e-7, δ=1.2)
    return δ * τ * (τ * tol/err)^r
end
