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
    println(io, tab(indent), Base.typename(typeof(prob)).wrapper, ":")
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
