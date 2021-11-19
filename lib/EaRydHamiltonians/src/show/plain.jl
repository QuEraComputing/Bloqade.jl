# Custom multi-line printing
_print(io::IO, x::AbstractFloat) = printstyled(io, round(x, sigdigits=3); color=:green)
_print(io::IO, x::Number) = printstyled(io, x; color=:green)

function _print(io::IO, x)
    if _iscallable(x)
        print(io, x, "(t)")
    else
        print(io, x)
    end
end

_print(io::IO, xs...) = foreach(x->_print(io, x), xs)

function _print_eachterm(f, io::IO, nsites::Int)
    indent = get(io, :indent, 0)
    limit = get(io, :limit, false)

    if limit && nsites > 6
        print(io, " "^indent);f(1);println(io, " +")
        print(io, " "^indent);f(2);println(io, " +")
        print(io, " "^indent);f(3);println(io, " +")

        println(io, " "^(indent + 2), "⋯")

        print(io, " "^indent);f(nsites-2);println(io, " +")
        print(io, " "^indent);f(nsites-1);println(io, " +")
        print(io, " "^indent);f(nsites)
    else
        for k in 1:nsites
            print(io, " "^indent)
            f(k)
            if k != nsites
                println(io, " +")
            end
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", t::AbstractTerm)
    indent = get(io, :indent, 0)
    println(io, " "^indent, nameof(typeof(t)))
    print_term(IOContext(io, :indent=>indent + 1), t)
end

print_term(io::IO, t::ZTerm) = _print_zterm(io, t.nsites, t.Δs)
print_term(io::IO, t::NTerm) = _print_nterm(io, t.nsites, t.Δs)
print_term(io::IO, t::XTerm) = _print_xterm(io, t.nsites, t.Ωs, t.ϕs)

function print_term(io::IO, t::RydInteract)
    indent = get(io, :indent, 0)
    _print_sum(io, nsites(t))
    _print(io, t.C)
    print(io, "/|r_i - r_j|^6 ")
    printstyled(io, "n_i n_j", color=:light_blue)
end

function print_term(io::IO, t::Hamiltonian)
    indent = get(io, :indent, 2)
    for (i, each) in enumerate(t.terms)
        println(io, " "^indent, " Term ", i)
        print_term(IOContext(io, :indent=>indent+2), each)
        
        if i != lastindex(t.terms)
            println(io)
            println(io)
        end
    end
end

function print_term(io::IO, t::Negative)
    print_term(IOContext(io, :negative=>true), t.term)
end

# NOTE: This should not be used in performance required code
# calculation intensive code should use generated version.
_iscallable(f) = !isempty(methods(f))

function _print_single_xterm(io::IO, Ω, ϕ)
    if _iscallable(Ω) || !(Ω ≈ 2)
        _print(io, Ω, "/2")
    end

    if !_iscallable(ϕ) && (isnothing(ϕ) || iszero(ϕ))
        printstyled(io, " σ^x", color=:light_blue)
    else
        _print(io, " (e^{", ϕ, "i}")
        printstyled(io, "|0)⟨1|", color=:light_blue)
        _print(io, " + e^{-", ϕ, "i}")
        printstyled(io, "|1⟩⟨0|", color=:light_blue)
        print(io, ")")
    end
end

function _print_sum(io::IO, nsites::Int)
    indent = get(io, :indent, 0)
    negative = get(io, :negative, false)
    print(io, " "^indent)
    negative && print(io, "-")
    print(io, "∑(n=1:$nsites) ")
end

function _print_xterm(io::IO, nsites::Int, Ω, ϕ)
    _print_sum(io, nsites)
    _print_single_xterm(io, Ω, ϕ)
end

function _print_xterm(io::IO, nsites::Int, Ω, ϕs::ConstParamListType)
    _print_eachterm(io, nsites) do k
        _print_single_xterm(io, Ω, getscalarmaybe(ϕs, k))
    end
end

function _print_xterm(io::IO, nsites::Int, Ωs::ConstParamListType, ϕ)
    _print_eachterm(io, nsites) do k
        _print_single_xterm(io, getscalarmaybe(Ωs, k), ϕ)
    end
end

function _print_xterm(io::IO, nsites::Int, Ωs::ConstParamListType, ϕs::ConstParamListType)
    _print_eachterm(io, nsites) do k
        _print_single_xterm(io, getscalarmaybe(Ωs, k), getscalarmaybe(ϕs, k))
    end
end

function _print_zterm(io::IO, nsites::Int, Δ)
    _print_sum(io, nsites)
    _print(io, Δ)
    printstyled(io, " σ^z", color=:light_blue)
end

function _print_zterm(io::IO, nsites::Int, Δs::ConstParamListType)
    _print_eachterm(io, nsites) do k
        _print(io, getscalarmaybe(Δs, k))
        printstyled(io, " σ^z", color=:light_blue)
    end
end

function _print_nterm(io::IO, nsites::Int, Δ)
    _print_sum(io, nsites)
    _print(io, Δ)
    printstyled(io, " n", color=:light_blue)
end

function _print_nterm(io::IO, nsites::Int, Δs::ConstParamListType)
    _print_eachterm(io, nsites) do k
        _print(io, getscalarmaybe(Δs, k))
        printstyled(io, " n", color=:light_blue)
    end
end
