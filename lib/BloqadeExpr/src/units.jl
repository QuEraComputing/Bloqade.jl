default_unit(unit, x) = x

default_unit(unit, x::Quantity) = uconvert(unit, x).val
default_unit(unit::typeof(NoUnits), x::Quantity) = uconvert(unit, x)
function default_unit(unit, xs::Tuple)
    return map(xs) do x
        default_unit(unit, x)
    end
end

function default_unit(unit, range::AbstractRange)
    a = default_unit(unit, first(range))
    b = default_unit(unit, step(range))
    c = default_unit(unit, last(range))
    return a:b:c
end

function default_unit(unit, x::AbstractArray{S}) where {T, S <: Quantity{T}}
    y = similar(x, T)
    @inbounds for i in eachindex(x)
        y[i] = default_unit(unit, x[i])
    end
    return y
end
