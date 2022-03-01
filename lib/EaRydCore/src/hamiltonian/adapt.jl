function Adapt.adapt_structure(to, t::XTerm)
    XTerm(t.nsites, adapt(to, t.Ωs), adapt(to, t.ϕs))
end

function Adapt.adapt_structure(to, t::ZTerm)
    ZTerm(t.nsites, adapt(to, t.Δs))
end

function Adapt.adapt_structure(to, t::NTerm)
    NTerm(t.nsites, adapt(to, t.Δs))
end

function Adapt.adapt_structure(to, t::RydInteract)
    RydInteract(adapt(to, t.atoms), adapt(to, t.C))
end

function Adapt.adapt_structure(to, t::Negative)
    Negative(adapt(to, t.term))
end

function Adapt.adapt_structure(to, t::Hamiltonian)
    Hamiltonian(map(x->Adapt.adapt(to, x), t.terms))
end
