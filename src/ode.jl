function Adapt.adapt_structure(to::CUDA.CuArrayAdaptor, cache::EquationCache{<:SparseMatrixCSC})
    # NOTE: cuSPARSE doesn't support complex-real mv routine
    # TODO: use a native Julia implementation?
    T = Complex{real(eltype(cache.hamiltonian))}
    return ContinuousEmulator.EquationCache(
        CuSparseMatrixCSR{T}(cache.hamiltonian),
        CuVector(cache.state),
    )
end

function RydbergEmulator.emulate!(r::Yao.ArrayReg{B, ST}, t::Real, h::AbstractTerm;
        algo=Vern8(), reltol=1e-8, abstol=1e-8, kwargs...
    ) where {B, ST <: CuArray}

    # NOTE: we force the hamiltonian matrix to be Complex on GPU
    # since cuSPARSE does not support complex-real mv routine
    precision_t = real(eltype(r.state))
    T = Complex{precision_t}
    t = precision_t(t)

    prob = ODEProblem(
        cu(ShordingerEquation(T, h)),
        vec(r.state), (zero(t), t);
        save_everystep=false, save_start=false, alias_u0=true, kwargs...
    )

    result = solve(prob, algo; reltol, abstol)
    return r
end

function RydbergEmulator.emulate!(r::RydbergReg{N, B, ST}, t::Real, h::AbstractTerm;
        algo=Vern8(), reltol=1e-8, abstol=1e-8, kwargs...
    ) where {N, B, ST <: CuArray}

    # NOTE: we force the hamiltonian matrix to be Complex on GPU
    # since cuSPARSE does not support complex-real mv routine
    precision_t = real(eltype(r.state))
    T = Complex{precision_t}
    t = precision_t(t)

    # convert subspace back to CPU, since we don't have a CUDA kernel for to_matrix!
    # and the performance of this function doesn't matter
    subspace = cpu(r.subspace)

    prob = ODEProblem(
        cu(ShordingerEquation(T, subspace, h)),
        vec(r.state), (zero(t), t);
        save_everystep=false, save_start=false, alias_u0=true, kwargs...
    )
    result = solve(prob, algo; reltol, abstol)
    return r
end
