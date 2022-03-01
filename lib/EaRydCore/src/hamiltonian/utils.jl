# TODO: use an actual symbolic zero, maybe consider SymbolicUtils Term
# struct Zero <: Number end

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
        y[i] = default_unit(unit, x[i]).val
    end
    return y
end

# function default_unit(unit, f) # parameters are function
#     return function waveform(t)
#         default_unit(unit, f(t))
#     end
# end

assert_has_time_method(::ConstParamType, name) = nothing
assert_has_time_method(::Nothing, name) = nothing # skip symbolic zero (currently nothing)
function assert_has_time_method(fs::Union{AbstractVector, Tuple}, name)
    for f in fs
        assert_has_time_method(f, name)
    end
end

function assert_has_time_method(f, name)
    hasmethod(f, Tuple{Real}) || throw(ArgumentError("invalid input for $name: method $f(::Real) is not defined"))
end

function assert_nsites(nsites::Int, p, name)
    p isa AbstractVector || p isa Tuple || return
    nsites == length(p) ||
        throw(ArgumentError(
            "nsites does not match size of $name " *
            "expect $nsites, got $(length(p))"
    ))
    return
end