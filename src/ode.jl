function Adapt.adapt_structure(to::CUDA.CuArrayAdaptor, cache::EquationCache{<:SparseMatrixCSC})
    # NOTE: cuSPARSE doesn't support complex-real mv routine
    # TODO: use a native Julia implementation?
    T = Complex{real(eltype(cache.hamiltonian))}
    return ContinuousEmulator.EquationCache(
        CuSparseMatrixCSR{T}(cache.hamiltonian),
        CuVector(cache.state),
    )
end

const CuRegister = Union{Yao.ArrayReg{B, ST}, RydbergReg{N, B, ST}} where {N, B, ST <: CuArray}

function ContinuousEmulator.ContinuousEvolution{P}(
    r::CuRegister, (start, stop)::Tuple{<:Real, <:Real}, h::AbstractTerm; kw...) where {P <: AbstractFloat}

    options = ContinuousOptions(;kw...)
    start = P(RydbergEmulator.default_unit(μs, start))
    stop = P(RydbergEmulator.default_unit(μs, stop))
    time = (start, stop)
    reg = adapt(RydbergEmulator.PrecisionAdaptor(P), r)
    space = RydbergEmulator.get_space(r)

    H = SparseMatrixCSC{Complex{P}}(h(start+1e-5), space)
    # NOTE: hamiltonian is Hermitian, so the nonzero position
    # is the same as transpose and we don't care about the value
    # of each entries for pre-allocate cache, they will be updated
    # to the right value later.
    dH = CuSparseMatrixCSR(
        CuVector{Cint}(H.colptr),
        CuVector{Cint}(rowvals(H)),
        CuVector{ComplexF32}(nonzeros(H)),
        size(H)
    )
    dstate = CuVector{Complex{P}}(undef, size(H, 1))
    dcache = ContinuousEmulator.EquationCache(dH, dstate)
    eq = ShordingerEquation(cu(space), h, dcache)

    ode_prob = ODEProblem(
        eq, vec(Yao.state(reg)), time;
        save_everystep=false, save_start=false, alias_u0=true,
        progress=options.progress,
        progress_steps=options.progress_steps,
    )
    return ContinuousEvolution{P}(reg, time, eq, ode_prob, options)
end
