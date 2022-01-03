"""
Built-in kernel functions.
"""
module Kernels

export gaussian,
    triangle,
    uniform,
    parabolic,
    biweight,
    tricube,
    cosine,
    logistic,
    sigmoid

# kernels
# https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use
_in_radius(t, value) = abs(t) ≤ 1 ? value : zero(t)
gaussian(t) = exp(-abs2(t)/2)/sqrt(2π)
triangle(t) = _in_radius(t, 1 - abs(t))

function uniform(t::T) where T
    # we calculate the normalize factor
    # in smooth, no need to div by m
    return T(abs(t) ≤ 1)
end

parabolic(t) = _in_radius(t, 3/4 * (1 - abs2(t)))
biweight(t) = _in_radius(t, 15/16 * (1 - abs2(t))^2)
triweight(t) = _in_radius(t, 35/32 * (1 - abs2(t))^3)
tricube(t) = _in_radius(t, 70/81 * (1 - abs(t)^3)^3)
cosine(t) = _in_radius(t, π/4 * cos(π/2 * t))
# TODO: check numerical stability
logistic(t) = 1 / (exp(t) + 2 + exp(-t))
sigmoid(t) = 2/(π * (exp(t) + exp(-t)))

end


# based on JuliaImages/edgepad
function edge_pad(vec::AbstractVector{T}, n::Int) where T
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
- `f`: a [`PiecewiseLinear`](@ref) function or a [`Waveform{<:PiecewiseLinear}`](@ref).

# Keyword Arguments

- `kernel_radius`: radius of the kernel.
- `edge_pad_size`: the size of edge padding.
"""
function smooth end

# forward waveform objects
function smooth(kernel, wf::Waveform{<:PiecewiseLinear}; kernel_radius::Real=0.3, edge_pad_size::Int=length(wf.f.clocks))
    return Waveform(smooth(kernel, wf.f; edge_pad_size, kernel_radius), wf.duration)
end

function smooth(wf::Waveform{<:PiecewiseLinear}; kernel_radius::Real=0.3, edge_pad_size::Int=length(wf.f.clocks))
    return smooth(Kernels.gaussian, wf; edge_pad_size, kernel_radius)
end

function smooth(kernel, f::PiecewiseLinear; kernel_radius::Real=0.3, edge_pad_size::Int=length(f.clocks))
    clocks = edge_pad(f.clocks, edge_pad_size)
    values = edge_pad(f.values, edge_pad_size)
    return smooth(kernel, clocks, values, kernel_radius)
end

smooth(f::PiecewiseLinear; kernel_radius::Real=0.3, edge_pad_size::Int=length(f.clocks)) =
    smooth(Kernels.gaussian, f; edge_pad_size, kernel_radius)

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
            !isequal(x, "eval")
        end |> xs->map(x->string("[`", x, "`](@ref)"), xs), "; "
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
            kernel(norm(x-x_i, 2)/kernel_radius)
        end

        Ŷ = sum(zip(Xi, Yi)) do (x_i, y_i)
            kernel(norm(x-x_i, 2)/kernel_radius) * y_i
        end
        return Ŷ / A
    end
end
