function Base.:(==)(x::RydInteract, y::RydInteract)
    return (x.atoms == y.atoms) && (x.C == y.C)
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

space_size(term::AbstractTerm, s::FullSpace) = 1 << nsites(term)
space_size(::AbstractTerm, s::Subspace) = length(s)