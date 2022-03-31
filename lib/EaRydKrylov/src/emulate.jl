@option struct KrylovOptions
    progress::Bool = false
    progress_step::Int = 1
    progress_name::String = "emulating"
    normalize_step::Int = 5
    normalize_finally::Bool = true
    tol::Float64 = 1e-7
end

struct KrylovEvolution{S, T <: Real, H <: Hamiltonian}
    reg::S
    start_clock::T
    durations::Vector{T}
    hamiltonian::H
    options::KrylovOptions
end

function KrylovEvolution(reg, h; start_clock=zero(eltype(durations)), kw...)
    options = from_kwargs(KrylovOptions; kw...)
    P = eltype(durations)
    T = isreal(h) ? P : Complex{P}
    return KrylovEvolution(reg, start_clock, durations, Hamiltonian(T, h), options)
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

function emulate!(prob::KrylovEvolution)
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
