Base.@propagate_inbounds function Base.iterate(prob::KrylovEvolution, step::Int=1)
    step > length(prob.durations) && return
    emulate_step!(prob, step)
    info = (step, prob.reg, duration=prob.durations[step], hamiltonian=prob.hs[step])
    return info, step+1
end

"""
    emulate_step!(prob::KrylovEvolution, step::Int=1)

Run single step evolution at given `step` index for a discrete evolution
(Krylov-based).
"""
Base.@propagate_inbounds function emulate_step!(prob::KrylovEvolution, step::Int=1)
    reg = prob.reg
    duration = prob.durations[step]
    h = prob.hs[step]
    emulate_routine!(reg, duration, h, prob.cache.H)

    if mod(step, prob.options.normalize_step) == 0
        normalize!(prob.reg)
    end

    if prob.options.normalize_finally && step == length(prob.durations)
        normalize!(prob.reg)
    end
    return prob
end

"""
    emulate!(evolution)

Run emulation on given evolution object.
See also [`KrylovEvolution`](@ref), or [`ContinousEvolution`](@ref).

# Arguments

- `evolution`: the evolution object.
"""
function emulate! end

function emulate!(prob::KrylovEvolution)
    niterations = length(prob.durations)
    @inbounds if prob.options.progress
        ProgressLogging.progress() do id
            for info in prob
                if prob.options.progress && mod(idx, prob.options.progress_step) == 0
                    @info prob.options.progress_name progress=info.step/niterations _id=id
                end
            end
        end
    else
        for info in prob; end
    end
    return prob
end

"""
    trotterize([start::Real=0], stop::Real, h::AbstractTerm; nsteps::Int=1000)

Trotterize time evolution of hamiltonian `h`. Each trotterize step
uses the first value of the interval as the clock to get `h(t)`, e.g
when start time is `0.1` and stop time is `0.5` with `1000` steps,
the first value of `h` will be `h(0.1)`, then `h(0.1 + dt)` and so on.

# Arguments

- `start`: start time, a real number, default is zero.
- `stop`: stop time, a real number.
- `h`: the hamiltonian expression.

# Keyword Arguments

- `nsteps`: number of steps in the evolution.
"""
function trotterize(start::Real, stop::Real, h::AbstractTerm; nsteps::Int=1000)
    schedule = range(;start, stop, length=nsteps)
    iteration = Iterators.take(schedule, nsteps-1)
    durations = [step(schedule) for _ in iteration]
    hs = [h(t) for t in iteration]
    return durations, hs
end

trotterize(duration::Real, h::AbstractTerm; nsteps::Int=1000) =
    trotterize(zero(duration), duration, h; nsteps)
