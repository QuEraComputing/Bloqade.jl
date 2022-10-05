assert_has_time_method(fs::Number, name) = nothing # skip scalars

function assert_has_time_method(fs::Vector, name)
    for f in fs
        assert_has_time_method(f, name)
    end
end

function assert_has_time_method(f, name)
    return is_time_function(f) || throw(ArgumentError("invalid input for $name: method $f(::Real) is not defined"))
end

function assert_nsites(nsites::Int, p, name)
    p isa AbstractVector || p isa Tuple || return
    nsites == length(p) ||
        throw(ArgumentError("nsites does not match size of $name " * "expect $nsites, got $(length(p))"))
    return
end

function assert_param(nsites::Int, p, name)
    assert_nsites(nsites, p, name)
    return assert_has_time_method(p, name)
end

# NOTE: is_time_function is only used in constructor
# don't use it in other places since it's slow

"""
    is_time_function(f)

Check if a function has method whose
signature is of `Tuple{Real}`
or `Tuple{T} where {T <: Real}`.
"""
function is_time_function(f)
    hasmethod(f, Tuple{Real}) && return true
    for T in concrete_subtypes(Real)
        hasmethod(f, Tuple{T}) && return true
    end
    return false
end

function is_const_param(x)
    return x isa Number || x isa Vector{<:Number}
end

concrete_subtypes(::Type{T}) where {T} = concrete_subtypes!([], T)

function concrete_subtypes!(list::Vector{Any}, ::Type{T}) where {T}
    isconcretetype(T) && return push!(list, T)
    for S in subtypes(T)
        isconcretetype(S) && push!(list, S)
        concrete_subtypes!(list, S)
    end
    return list
end

is_time_function(fs::Vector) = all(is_time_function, fs)
