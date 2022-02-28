const CuRegister = Union{Yao.ArrayReg{1,T, State},RydbergReg{Layout,State}} where {T, Layout,State<:CuArray}

function EaRydODE.ODEEvolution{P}(
    r::CuRegister, (start, stop)::Tuple{<:Real,<:Real}, h::AbstractTerm; kw...
) where {P<:AbstractFloat}
    layout = EaRydCore.MemoryLayout(r)
    if layout isa RealLayout
        isreal(h) || error("cannot use RealLayout for non-real hamiltonian")
    end

    main_options, ode_options = [], []
    for (k, v) in kw
        if k in fieldnames(ODEOptions)
            push!(main_options, k=>v)
        else
            push!(ode_options, k=>v)
        end
    end

    options = ODEOptions(;main_options...)
    start = EaRydCore.default_unit(μs, start)
    stop = EaRydCore.default_unit(μs, stop)
    time = (start, stop)
    reg = adapt(EaRydCore.PrecisionAdaptor(P), r)
    space = EaRydCore.get_space(r)

    if layout isa RealLayout && isreal(h)
        T = P
    else
        T = Complex{P}
    end

    H = SparseMatrixCSC{T,Cint}(h(start + sqrt(eps(P))), space)
    # NOTE: hamiltonian is Hermitian, so the nonzero position
    # is the same as transpose and we don't care about the value
    # of each entries for pre-allocate cache, they will be updated
    # to the right value later.
    dH = CuSparseMatrixCSR(
        CuVector{Cint}(H.colptr), CuVector{Cint}(rowvals(H)), CuVector{T}(nonzeros(H)), size(H)
    )

    if layout isa ComplexLayout
        dstate = CuVector{Complex{P}}(undef, size(H, 1))
    elseif layout isa RealLayout
        dstate = CuMatrix{P}(undef, size(H, 1), 2)
    else
        error("$layout is not supported")
    end

    dcache = EaRydODE.EquationCache(dH, layout, dstate)
    eq = SchrodingerEquation(h, cu(space), dcache)

    ode_prob = ODEProblem(
        eq,
        Yao.statevec(reg),
        time;
        save_everystep=false,
        save_start=false,
        alias_u0=true,
        progress=options.progress,
        progress_name=options.progress_name,
        progress_steps=options.progress_steps,
        ode_options...
    )
    return ODEEvolution{P}(reg, time, eq, ode_prob, options)
end
