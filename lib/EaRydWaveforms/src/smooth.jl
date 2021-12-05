"""
Built-in kernel functions.
"""
module Kernels

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

parabolic(t) = _in_radius(3/4 * (1 - abs2(t)))
biweight(t) = _in_radius(15/16 * (1 - abs2(t))^2)
triweight(t) = _in_radius(35/32 * (1 - abs2(t))^3)
tricube(t) = _in_radius(70/81 * (1 - abs(t)^3)^3)
cosine(t) = _in_radius(π/4 * cos(π/2 * t))
# TODO: check numerical stability
logistic(t) = 1 / (exp(t) + 2 + exp(-t))
sigmoid(t) = 2/(π * (exp(u) + exp(-u)))

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

# with gaussian filter
function smooth(filter, wf::Waveform{<:PiecewiseLinear})
    clocks, values = wf.f.clocks, wf.f.values
    return Waveform(
        smooth(filter, clocks, values, radius),
        wf.interval
    )
end

"""
    smooth(kernel, Xi::Vector, Yi::Vector, kernel_radius::Real)

Kernel smoother function.

!!! tip

    This method is slower when the number of data points
    become large because this calculates only a single
    point at each time. This is only useful when one does
    not know the input point `x`, consider to use
    [`fast_smooth`](@ref) when the desired input `X`
    sequence is known.

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

# # this is an optimized version of above using FFT
# function fast_smooth(kernel, X::Vector, Xi::Vector, Yi::Vector, kernel_radius::Real)
#     x_conv = [kernel(norm(x - x_i, 2)/kernel_radius) for (x, x_i) in zip(Xi, Yi)]
# end
