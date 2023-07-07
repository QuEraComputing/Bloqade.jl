"""
    struct KrylovOptions

Krylov evolution options.
"""
@option struct KrylovOptions
    progress::Bool = false
    progress_step::Int = 1
    progress_name::String = "emulating"
    normalize_step::Int = 5
    normalize_finally::Bool = true
    tol::Float64 = 1e-7
    expmv_backend = expmv!
end

## this is the abstract type for all evolvers 
## KrylovEvolution, Magnus4Evolution, CFET42Evolution
abstract type Evolver end


# These are common parts shared by all the evolvers
Base.length(prob::Evolver) = length(prob.durations) + 1

function Base.iterate(prob::Evolver)
    info = (; step = 1, reg = prob.reg, clock = prob.start_clock, duration = zero(prob.start_clock))
    return info, (2, prob.start_clock)
end

Base.@propagate_inbounds function Base.iterate(prob::Evolver, (step, clock))
    step > length(prob) && return

    duration = prob.durations[step-1]
    emulate_step!(prob, step, clock, duration)

    info = (; step, reg = prob.reg, clock = clock + duration, duration)
    return info, (step + 1, clock + duration)
end

## driver function, user entry point
function BloqadeExpr.emulate!(prob::Evolver)
    niterations = length(prob)
    @inbounds if prob.options.progress
        ProgressLogging.progress() do id
            for info in prob
                if prob.options.progress && mod(info.step, prob.options.progress_step) == 0
                    @info prob.options.progress_name progress = info.step / niterations _id = id
                end
            end
        end
    else
        for info in prob
        end
    end
    return prob
end

tab(indent) = " "^indent

function Base.show(io::IO, mime::MIME"text/plain", prob::Evolver)
    indent = get(io, :indent, 0)
    if typeof(prob) <: CFETEvolution
        println(io, tab(indent), "CFETEvolution", "<",Base.typename(typeof(prob.alg_table)).wrapper, ">:")
    else
        println(io, tab(indent), Base.typename(typeof(prob)).wrapper, ":")
    end
    # state info
    print_state_info(io, prob)
    println(io)

    # clocks
    println(io, tab(indent + 2), "clocks")
    println(io, tab(indent + 4), "start:", prob.start_clock, "μs")
    println(io, tab(indent + 4), " last:", prob.start_clock + sum(prob.durations), "μs")
    println(io)

    # equation info
    show(IOContext(io, :indent => indent + 2), mime, prob.hamiltonian)
    println(io)
    println(io)

    println(io, tab(indent + 2), "Options:")
    for name in fieldnames(KrylovOptions)
        println(io, tab(indent + 4), name, "=", repr(getfield(prob.options, name)))
    end
end

function print_state_info(io::IO, prob::Evolver)
    indent = get(io, :indent, 0)
    println(io, tab(indent + 2), "register info:")
    print(io, tab(indent + 4), "type: ")
    printstyled(io, typeof(prob.reg); color = :green)
    println(io)

    print(io, tab(indent + 4), "storage size: ")
    printstyled(io, Base.format_bytes(storage_size(prob.reg)); color = :yellow)
    return println(io)
end


function uniquetol_list(itr::Array{T,1},tol=1e-6) where {T<:Complex}
    preal = sortperm(real(itr))     # buffer for permutation idxs of reals
    dreal = diff(real(itr[preal]))  # buffer for differences of reals

    pimag = Array{Int}(size(preal))     # buffer for permutation idxs of imags
    dimag = Array{Float64}(size(dreal)) # buffer for differences of imags

    out = [ ];
    pri1=1;
    for pri2 = 1:length(preal)
        # look for tol-wide gap in reals
        if pri2==length(preal) || dreal[pri2] > tol
            pridxs = preal[pri1:pri2]
            sortperm!( pimag[1:length(pridxs)], imag(itr[pridxs]) )
            dimag[1:(length(pridxs)-1)] = diff(imag(itr[preal[pridxs]]))

            pii1=1;
            for pii2 = 1:length(pridxs)
                if pii2==length(pridxs) || dimag[pii2] > tol
                    push!(out, Array{T,1}(itr[pridxs[pii1:pii2]]))
                    pii1=pii2+1;
                end # imag tol check
            end # pii2, tol-large gap in imags

            pri1=pri2+1;
        end # real tol check
    end # pri2, tol-large gap in reals
    return out
end