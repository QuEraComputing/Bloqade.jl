using Printf

abstract type AbstractTerm end

struct InteractTerm{T, Atom <: RydAtom} <: AbstractTerm
    C::T
    atoms::Vector{Atom}
end

struct XTerm{Omega, Phi} <: AbstractTerm
    Ωs::Omega
    ϕs::Phi
end

struct ZTerm{Delta} <: AbstractTerm
    Δs::Delta
end

struct Hamiltonian{Terms <: Tuple} <: AbstractTerm
    terms::Terms
end

XTerm(Ωs) = XTerm(Ωs, nothing)

function Base.show(io::IO, t::InteractTerm)
    print(io, "C/|r_i - r_j|^6 ")
    printstyled(io, "n_i n_j", color=:light_blue)
end

function Base.show(io::IO, t::XTerm)
    print(io, "Ω ⋅ (e^{iϕ}")
    printstyled(io, "|0)", color=:light_blue)
    print(io, "+e^{-iϕ}")
    printstyled(io, "|1⟩", color=:light_blue)
end

function Base.show(io::IO, t::XTerm{<:Any, Nothing})
    print(io, "Ω")
    printstyled(io, " σ^x", color=:light_blue)
end

function Base.show(io::IO, t::ZTerm)
    print(io, "Δ")
    printstyled(io, " σ^z", color=:light_blue)
end

function Base.show(io::IO, x::Hamiltonian)
    print(io, x.terms[1])
    for t in x.terms[2:end]
        print(io, " + ")
        print(io, t)
    end
end

function Base.:(+)(x::AbstractTerm, y::AbstractTerm)
    return Hamiltonian((x, y))
end

function Base.:(+)(x::AbstractTerm, y::Hamiltonian)
    return Hamiltonian((x, y.terms...))
end

function Base.:(+)(x::Hamiltonian, y::AbstractTerm)
    return Hamiltonian((x.terms..., y))
end

function update!(t::XTerm, Ωs, ϕs) 
end

function update!(t::ZTerm, Δs)
end

function update!(t::InteractTerm, C, atoms)
end

h = InteractTerm(2.2, rand(10)) + ZTerm(rand(10)) + XTerm(rand(10))
h = InteractTerm(2.2, rand(10)) + ZTerm(rand(10)) + XTerm(rand(10), rand(10))
