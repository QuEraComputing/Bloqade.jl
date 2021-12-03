struct CompositeWaveform{N, T, Waveforms <: Tuple} <: AbstractWaveform
    waveforms::Waveforms
    durations::NTuple{N, T}
    checkpoints::NTuple{N, T}
end

"""
    CompositeWaveform(waveforms::AbstractWaveform...)

Composite waveforms into a new sequence.

# Example

```julia
julia> waveform = CompositeWaveform(
               SinusoidalWaveform(duration=2.2),
               RampWaveform(;duration=0.5, start=0.0, stop=1.0)
           )
                  ⠀⠀⠀⠀⠀⠀⠀⠀CompositeWaveform{Float64}⠀⠀⠀⠀⠀⠀⠀⠀ 
                  ┌────────────────────────────────────────┐ 
                1 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠊⠉⠉⠉⠉⠑⠢⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠓⡄⠀⠀⠀⠀⢰⠁⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡰⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⡜⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⢀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⡜⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⢀⠇⠀⠀⠀⠀⠀⠀│ 
   value (rad/µs) │⠀⠀⠀⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⡸⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⢰⠁⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⢀⡎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣧⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                0 │⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠟⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  └────────────────────────────────────────┘ 
                  ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀3⠀ 

```
"""
function CompositeWaveform(waveforms::AbstractWaveform...)
    durations = map(duration, waveforms)
    return CompositeWaveform(waveforms, durations, accumulate(+, durations))
end

Base.eltype(waveform::CompositeWaveform) = eltype(waveform.waveforms[1])
Base.getindex(waveform::CompositeWaveform, idx::Int) = waveform.waveforms[idx]
Base.length(waveform::CompositeWaveform) = length(waveform.waveforms)

duration(waveform::CompositeWaveform) = waveform.checkpoints[end]

function sample_values(waveform::CompositeWaveform, dt::Real)
    values = typeof(dt)[]
    for wf in waveform.waveforms
        for clock in 0:dt:duration(wf)-dt
            push!(values, wf(clock))
        end
    end

    # push last point to match the inteval
    # length
    wf = last(waveform.waveforms)
    push!(values, wf(duration(wf)))
    return values
end

# NOTE: change this to a generated function if the performance matters
# in the future
function (waveform::CompositeWaveform)(t::Real, offset::Real=zero(t))
    clock = t - offset
    idx = 1
    while idx < length(waveform) && clock > waveform.checkpoints[idx]
        idx += 1
    end
    offset = isone(idx) ? offset : offset + waveform.checkpoints[idx-1]
    return waveform.waveforms[idx](t, offset)
end
