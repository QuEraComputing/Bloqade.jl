@option struct KrylovOptions
    progress::Bool = false
    progress_step::Int = 1
    progress_name::String = "emulating"
    normalize_step::Int = 5
    normalize_finally::Bool = true
    tol::Float64 = 1e-7
end

struct KrylovEvolution{Reg <: AbstractRegister, T <: Real, H <: Hamiltonian}
    reg::Reg
    start_clock::T
    durations::Vector{T}
    hamiltonian::H
    options::KrylovOptions

    function KrylovEvolution{Reg, T, H}(reg, start_clock, durations, hamiltonian, options) where {Reg, T, H}
        start_clock ≥ 0 || throw(ArgumentError("start clock must not be negative"))
        all(≥(0), durations) || throw(ArgumentError("durations must not be negative"))
        new{Reg, T, H}(reg, start_clock, durations, hamiltonian, options)
    end
end

function KrylovEvolution(reg, start_clock, durations, hamiltonian, options)
    return KrylovEvolution{typeof(reg), typeof(start_clock), typeof(hamiltonian)}(
        reg, start_clock, durations, hamiltonian, options
    )
end

function Adapt.adapt_structure(to, x::KrylovEvolution)
    return KrylovEvolution(
        adapt(to, x.reg),
        x.start_clock, x.durations,
        adapt(to, x.hamiltonian),
        x.options
    )
end

function KrylovEvolution(reg::AbstractRegister, clocks, h; kw...)
    all(≥(0), clocks) || throw(ArgumentError("clocks must not be negative"))
    options = from_kwargs(KrylovOptions; kw...)
    P = real(eltype(statevec(reg)))
    T = isreal(h) ? P : Complex{P}
    start_clock, durations = first(clocks), diff(clocks)
    return KrylovEvolution(reg, start_clock, durations, Hamiltonian(T, h, space(reg)), options)
end

function emulate_step!(prob::KrylovEvolution, step::Int, clock::Real, duration::Real)
    state = statevec(prob.reg)
    h = prob.hamiltonian

    expmv!(-duration * im, h(clock), state; prob.options.tol)
    if mod(step, prob.options.normalize_step) == 0
        normalize!(prob.reg)
    end

    if prob.options.normalize_finally && step == length(prob.durations)
        normalize!(prob.reg)
    end
    return prob
end

Base.@propagate_inbounds function Base.iterate(prob::KrylovEvolution, (step, clock)=(1, prob.start_clock))
    step > length(prob.durations) && return

    duration=prob.durations[step]
    emulate_step!(prob, step, clock, duration)

    info = (;step, reg=prob.reg, duration)
    return info, (step+1, clock+duration)
end

function BloqadeExpr.emulate!(prob::KrylovEvolution)
    niterations = length(prob.durations)
    @inbounds if prob.options.progress
        ProgressLogging.progress() do id
            for info in prob
                if prob.options.progress && mod(info.step, prob.options.progress_step) == 0
                    @info prob.options.progress_name progress=info.step/niterations _id=id
                end
            end
        end
    else
        for info in prob; end
    end
    return prob
end
