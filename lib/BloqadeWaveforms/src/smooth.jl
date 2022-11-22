"""
Built-in kernel functions.
"""
module Kernels

export gaussian, triangle, uniform, parabolic, biweight, triweight, tricube, cosine, logistic, sigmoid

# kernels
# https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use
_in_radius(t, value) = abs(t) ≤ 1 ? value : zero(t)
"""
    gaussian(t)

Gaussian kernel function for smoothing waveforms via [`smooth`](@ref).

The function is defined as:
```math
f(t) = \\frac{1}{\\sqrt{2π}}e^{-\\frac{1}{2}|t|^2}
```
"""
gaussian(t) = exp(-abs2(t) / 2) / sqrt(2π)
"""
    triangle(t)

Triangle kernel function for smoothing waveforms via [`smooth`](@ref).

The function is defined as:
```math
f(t) = 1 - |t|
```
where ``|t| ≤ 1``. Otherwise, ``f(t) = 0``
"""
triangle(t) = _in_radius(t, 1 - abs(t))

"""
    uniform(t::T) where {T}

Uniform kernel function for smoothing waveforms via [`smooth`](@ref)

The function returns ``1`` for ``|t| ≤ 1``. Otherwise, it returns ``0``.
"""
function uniform(t::T) where {T}
    # we calculate the normalize factor
    # in smooth, no need to div by m
    return T(abs(t) ≤ 1)
end

"""
    parabolic(t)

Parabolic kernel function for smoothing waveforms via [`smooth`](@ref)

The function is defined as:
```math
f(t) = \\frac{3}{4}(1 - |t|^2)
```
when ``|t| ≤ 1``. Otherwise, ``f(t) = 0``.
"""
parabolic(t) = _in_radius(t, 3 / 4 * (1 - abs2(t)))
"""
    biweight(t)

Biweight kernel function for smoothing waveforms via [`smooth`](@ref)

The function is defined as:
```math 
f(t) = \\frac{15}{16}(1 - |t|^2)^2
```
when ``|t| ≤ 1``. Otherwise, ``f(t) = 0``.
"""
biweight(t) = _in_radius(t, 15 / 16 * (1 - abs2(t))^2)
"""
    triweight(t)

Biweight kernel function for smoothing waveforms via [`smooth`](@ref)

The function is defined as:
```math
f(t) = \\frac{35}{32}(1 - |t|^2)^3
```
when ``|t| ≤ 1``. Otherwise, ``f(t) = 0``.
"""
triweight(t) = _in_radius(t, 35 / 32 * (1 - abs2(t))^3)
"""
    tricube(t)

Tricube kernel function for smoothing waveforms via [`smooth`](@ref)

The function is defined as:
```math
f(t) = \\frac{70}{81}(1 - |t|^3)^3
```
when ``|t| ≤ 1``. Otherwise, ``f(t) = 0``.
"""
tricube(t) = _in_radius(t, 70 / 81 * (1 - abs(t)^3)^3)
"""
    cosine(t)

Cosine kernel function for smoothing waveforms via [`smooth`](@ref)

The function is defined as:
```math
f(t) = \\frac{π}{4}\\cos\\left(\\frac{π}{2}t\\right)
```
when ``|t| ≤ 1``. Otherwise, ``f(t) = 0``.
"""
cosine(t) = _in_radius(t, π / 4 * cos(π / 2 * t))
# TODO: check numerical stability
"""
    logistic(t)

Logistic kernel function for smoothing waveforms via [`smooth`](@ref)

The function is defined as:
```math
f(t) = \\frac{1}{e^{t} + 2 + e^{-t}}
```
"""
logistic(t) = 1 / (exp(t) + 2 + exp(-t))
"""
    simgoid(t)

Sigmoid kernel funciton for smoothing waveforms via [`smooth`](@ref)

The function is defined as:
```math
f(t) = \\frac{2}{\\pi (e^{t} + e^{-t})}
```
"""
sigmoid(t) = 2 / (π * (exp(t) + exp(-t)))

end

# based on JuliaImages/edgepad
function edge_pad(vec::AbstractVector{T}, n::Int) where {T}
    res = zeros(T, length(vec) + 2n)
    @inbounds begin
        for idx in 1:n
            res[idx] = vec[1]
            res[end-idx+1] = vec[end]
        end
        res[n+1:end-n] .= vec
    end
    return res
end

"""
    smooth([kernel=Kernel.gaussian], f; edge_pad_size::Int=length(f.clocks))

Kernel smoother function for piece-wise linear function/waveform via weighted moving average method.

# Arguments

- `kernel`: the kernel function, default is [`Kernels.gaussian`](@ref).
- `f`: a [`Union{PiecewiseLinear, PiecewiseConstant}`](@ref) function or a [`Waveform{<:Union{PiecewiseLinear, PiecewiseConstant}}`](@ref).

# Keyword Arguments

- `kernel_radius`: radius of the kernel.
- `edge_pad_size`: the size of edge padding.
"""
function smooth end

# forward waveform objects
function smooth(
    kernel,
    wf::Waveform{<:Union{PiecewiseLinear,PiecewiseConstant}};
    kernel_radius::Real = 0.3,
    edge_pad_size::Int = length(wf.f.clocks),
    dt = kernel_radius / 1e2,
)
    return Waveform(smooth(kernel, wf.f; edge_pad_size, kernel_radius, dt), wf.duration)
end

function smooth(
    wf::Waveform{<:Union{PiecewiseLinear,PiecewiseConstant}};
    kernel_radius::Real = 0.3,
    edge_pad_size::Int = length(wf.f.clocks),
    dt = kernel_radius / 1e2,
)
    return smooth(Kernels.gaussian, wf; edge_pad_size, kernel_radius, dt)
end

function smooth(
    kernel,
    f::Union{PiecewiseLinear,PiecewiseConstant};
    kernel_radius::Real = 0.3,
    edge_pad_size::Int = length(f.clocks),
    dt = kernel_radius / 1e2,
)
    clocks = zero(dt):dt:last(f.clocks)
    clocks = edge_pad(clocks, edge_pad_size)
    values = [f(t) for t in clocks]
    values = edge_pad(values, edge_pad_size)
    return smooth(kernel, clocks, values, kernel_radius)
end

smooth(
    f::Union{PiecewiseLinear,PiecewiseConstant};
    kernel_radius::Real = 0.3,
    edge_pad_size::Int = length(f.clocks),
    dt = kernel_radius / 1e2,
) = smooth(Kernels.gaussian, f; edge_pad_size, kernel_radius)

"""
    smooth(kernel, Xi::Vector, Yi::Vector, kernel_radius::Real)

Kernel smoother function via weighted moving average method.
See also [Kernel Smoother](https://en.wikipedia.org/wiki/Kernel_smoother).

# Theory

Kernel function smoothing is a technique to define
a smooth function ``f: \\mathcal{R}^p → \\mathbf{R}``
from a set of discrete points by weighted averaging
the neighboring points. It can be written as the following
equation.

```math
Ŷ(X) = \\sum_i K(X, X_i) Y_i / \\sum_i K(X, X_i)
```

where ``Ŷ(X)`` is the smooth function by calculating
the moving average of known data points ``X_i`` and ``Y_i``.
`K` is the kernel function, where ``K(\\frac{||X - X_i||}{h_λ})``
decrease when the Euclidean norm ``||X - X_i||`` increase,
``h_λ`` is a parameter controls the radius of the kernel.

# Available Kernels

The following kernel functions are available via the
[`Kernels`](@ref) module:

$(
    join(
        filter(names(Kernels; all=true)) do name
            x = string(name)
            !startswith(x, '_') &&
            !startswith(x, '#') &&
            !isequal(x, "Kernels") &&
            !isequal(x, "eval") &&
            !isequal(x, "include")
        end |> xs->map(x->string("[`Kernels.", x, "`](@ref)"), xs), "; "
    )
)

# Arguments

- `kernel`: a Julia function that has method `kernel(t::Real)`.
- `Xi::Vector`: a list of inputs `X_i`.
- `Yi::Vector`: a list of outputs `Y_i`.
- `kernel_radius::Real`: the radius of the kernel.
"""
function smooth(kernel, Xi::Vector, Yi::Vector, kernel_radius::Real)
    return function smoothed_function(x)
        A = sum(Xi) do x_i
            return kernel(norm(x - x_i, 2) / kernel_radius)
        end

        Ŷ = sum(zip(Xi, Yi)) do (x_i, y_i)
            return kernel(norm(x - x_i, 2) / kernel_radius) * y_i
        end
        return Ŷ / A
    end
end
