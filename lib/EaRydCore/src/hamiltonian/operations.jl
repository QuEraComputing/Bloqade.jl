# Base.:(+)(x::AbstractTerm, y::AbstractTerm) = Hamiltonian((x, y))
# Base.:(+)(x::AbstractTerm, y::Hamiltonian) = Hamiltonian((x, y.terms...))
# Base.:(+)(x::Hamiltonian, y::AbstractTerm) = Hamiltonian((x.terms..., y))
# Base.:(+)(x::Hamiltonian, y::Hamiltonian) = Hamiltonian((x.terms..., y.terms...))

# # absorb - to RHS
# Base.:(-)(x::AbstractTerm, y::AbstractTerm) = Hamiltonian((x, -y))
# Base.:(-)(x::AbstractTerm, y::Hamiltonian) = Hamiltonian((x, map(-, y.terms)...))
# Base.:(-)(x::Hamiltonian, y::AbstractTerm) = Hamiltonian((x.terms..., -y))
# Base.:(-)(x::Hamiltonian, y::Hamiltonian) = Hamiltonian((x.terms..., map(-, y.terms)...))

# Base.:(-)(x::XTerm{<:ConstParamType}) = XTerm(x.nsites, map(-, x.Ωs), x.ϕs)
# Base.:(-)(x::XTerm) = XTerm(x.nsites, t->-x.Ωs(t), x.ϕs)

# Base.:(-)(x::ZTerm{<:ConstParamType}) = ZTerm(x.nsites, map(-, x.Δs))
# Base.:(-)(x::ZTerm) = ZTerm(x.nsites, t->-x.Δs(t))
# Base.:(-)(x::NTerm{<:ConstParamType}) = NTerm(x.nsites, map(-, x.Δs))
# Base.:(-)(x::NTerm) = NTerm(x.nsites, t->-x.Δs(t))


function Base.:(==)(x::RydInteract, y::RydInteract)
    return (x.atoms == y.atoms) && (x.C == y.C)
end

function Base.:(==)(x::Hamiltonian, y::Hamiltonian)
    return all(t in y.terms for t in x.terms) && all(t in x.terms for t in y.terms)
end

"""
    getterm(terms, k, k_site)

Get the value of k-th local term in `terms`
given the site configuration as `k_site`.
"""
function getterm end

getterm(t::Negative, k, k_site) = -getterm(t.term, k, k_site)

function getterm(t::XTerm, k, k_site)
    if k_site == 0
        return getscalarmaybe(t.Ωs, k)/2 * exp(im * getscalarmaybe(t.ϕs, k))
    else
        return getscalarmaybe(t.Ωs, k)/2 * exp(-im * getscalarmaybe(t.ϕs, k))
    end
end

function getterm(t::XTerm{<:Any, Nothing}, k, k_site)
    return getscalarmaybe(t.Ωs, k)/2
end

function getterm(t::ZTerm, k, k_site)
    if k_site == 0
        return getscalarmaybe(t.Δs, k)
    else
        return -getscalarmaybe(t.Δs, k)
    end
end

function getterm(t::NTerm, k, k_site)
    if k_site == 1
        return getscalarmaybe(t.Δs, k)
    else
        return 0
    end
end

function getterm(t::Hamiltonian, k, k_site)
    error("composite Hamiltonian term cannot be indexed")
end

space_size(term::AbstractTerm, s::FullSpace) = 1 << nsites(term)
space_size(::AbstractTerm, s::Subspace) = length(s)