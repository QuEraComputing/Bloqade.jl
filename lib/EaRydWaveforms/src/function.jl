"""
    FunctionWaveform{F, T <: Real} <: AbstractWaveform
    
Type for waveform defined by a Julia function/callable object.
"""
struct FunctionWaveform{F, T <: Real} <: AbstractWaveform
    f::F
    duration::T
end

"""
    FunctionWaveform(f; duration::Real)

Create a `FunctionWaveform` from given function `f`
and `duration` of the waveform. See also [`SinusoidalWaveform`](@ref),
[`ConstantWaveform`](@ref).

# Example

```julia
julia> FunctionWaveform(sin, duration=4π)
                  ⠀⠀FunctionWaveform{typeof(sin), Float64}⠀⠀ 
                  ┌────────────────────────────────────────┐ 
                1 │⠀⠀⡞⠹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⠋⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⢸⠁⠀⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⡞⠀⠸⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⡞⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⢠⠇⠀⠀⣇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⢀⡇⠀⠀⠈⡇⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⢸⠀⠀⠀⠀⢧⠀⠀⠀⠀⠀⠀⠀⡏⠀⠀⠀⠘⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⡼⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⢀⡇⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⡇⠀⠀⠀⠀⠘⡆⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⢹⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   value (rad/µs) │⠧⠤⠤⠤⠤⠤⡧⠤⠤⠤⠤⠤⡼⠤⠤⠤⠤⠤⢼⠤⠤⠤⠤⠤⠤⡤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│ 
                  │⠀⠀⠀⠀⠀⠀⢳⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠈⡇⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⢠⠇⠀⠀⠀⠀⠀⠀⣇⠀⠀⠀⠀⣸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠈⡇⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⣇⠀⠀⠀⡞⠀⠀⠀⠀⠀⠀⠀⠸⡄⠀⠀⢠⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⢀⡇⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⣸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠈⡇⠀⣸⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
               -1 │⠀⠀⠀⠀⠀⠀⠀⠀⢹⣠⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣇⣸⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  └────────────────────────────────────────┘ 
                  ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀20⠀ 
```
"""
function FunctionWaveform(f; duration::Real)
    return FunctionWaveform(f, duration)
end

duration(wf::FunctionWaveform) = wf.duration
Base.eltype(::FunctionWaveform{<:Any, T}) where T = T

function (wf::FunctionWaveform)(t::Real, offset::Real=zero(t))
    assert_clock(t, duration(wf), offset)
    return wf.f(t - offset)
end

# this is for pretty printing & dispatch
struct Sinusoidal{T} <: Function
    amp::T
end

(f::Sinusoidal)(t::Real) = f.amp * sin(t)

struct Constant{T} <: Function
    amp::T
end

(f::Constant)(::Real) = f.amp

"""
    const SinusoidalWaveform{T} = FunctionWaveform{Sinusoidal{T}, T}
    SinusoidalWaveform(;duration::Real, amplitude::Real=1.0)

A convenient type alias for sinusoidal waveform. It is equivalent
to `FunctionWaveform(t->amplitude * sin(t), duration)`.

# Example

```julia
julia> SinusoidalWaveform(duration=4π, amplitude=2.1)
                  ⠀⠀⠀⠀⠀⠀⠀⠀SinusoidalWaveform{Float64}⠀⠀⠀⠀⠀⠀⠀ 
                  ┌────────────────────────────────────────┐ 
                3 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⡴⠲⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠖⢦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⣸⠁⠀⢳⠀⠀⠀⠀⠀⠀⠀⠀⠀⡏⠀⠘⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⢀⡇⠀⠀⠘⡆⠀⠀⠀⠀⠀⠀⠀⢸⠁⠀⠀⢹⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⢸⠀⠀⠀⠀⢳⠀⠀⠀⠀⠀⠀⠀⡏⠀⠀⠀⠘⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⡏⠀⠀⠀⠀⢸⡀⠀⠀⠀⠀⠀⢰⠃⠀⠀⠀⠀⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   value (rad/µs) │⠧⠤⠤⠤⠤⠤⡧⠤⠤⠤⠤⠤⡼⠤⠤⠤⠤⠤⢼⡤⠤⠤⠤⠤⢤⡤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│ 
                  │⠀⠀⠀⠀⠀⠀⢹⠀⠀⠀⠀⢀⡇⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠘⡆⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⢹⠀⠀⠀⠀⡏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⢣⠀⠀⠀⡏⠀⠀⠀⠀⠀⠀⠀⠘⡆⠀⠀⢰⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠘⡆⠀⢸⠁⠀⠀⠀⠀⠀⠀⠀⠀⢳⠀⠀⡞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠹⠤⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠧⠼⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
               -3 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  └────────────────────────────────────────┘ 
                  ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀20⠀ 
```
"""
const SinusoidalWaveform{T} = FunctionWaveform{Sinusoidal{T}, T}

"""
    const ConstantWaveform{T} = FunctionWaveform{Constant{T}, T}
    ConstantWaveform(;duration::Real, amplitude::Real=1.0)

A convenient type alias for constant waveform. It is equivalent
to `FunctionWaveform(t->amplitude, duration)`.

# Example

```julia
julia> ConstantWaveform(duration=4π, amplitude=1.2)
                  ⠀⠀⠀⠀⠀⠀⠀⠀⠀ConstantWaveform{Float64}⠀⠀⠀⠀⠀⠀⠀⠀ 
                  ┌────────────────────────────────────────┐ 
                3 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   value (rad/µs) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                0 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  └────────────────────────────────────────┘ 
                  ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀20⠀ 
```
"""
const ConstantWaveform{T} = FunctionWaveform{Constant{T}, T}

function SinusoidalWaveform(;duration::Real, amplitude::Real=1.0)
    return FunctionWaveform(Sinusoidal(amplitude), duration)
end

function ConstantWaveform(;duration::Real, amplitude::Real=1.0)
    return FunctionWaveform(Constant(amplitude), duration)
end
