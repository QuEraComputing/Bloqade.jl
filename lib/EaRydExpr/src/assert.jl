assert_has_time_method(fs::Number, name) = nothing # skip scalars

function assert_has_time_method(fs::Vector, name)
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

function assert_param(nsites::Int, p, name)
    assert_nsites(nsites, p, name)
    assert_has_time_method(p, name)
end

is_time_function(f) = hasmethod(f, Tuple{Real})
is_time_function(fs::Vector) = all(is_time_function, fs)
