# a linear map for low-level hamiltonian representation
struct Hamiltonian{FS <: Tuple, TS <: Tuple}
    fs::FS # prefactor of each term
    ts::TS # const linear map of each term
end

struct StepHamiltonian{T, FS, TS}
    t::T # clock
    h::Hamiltonian{FS, TS}
end

(h::Hamiltonian)(t::Real) = StepHamiltonian(t, h)

# NOTE: these high level expression
# will be converted to lower level Hamiltonian object
# before simulation no need to aggressively specialize on them
# the following operators are mainly used to interface
# with Yao

struct XPhase{T} <: PrimitiveBlock{2}
    ϕ::T
end

abstract type AbstractTerm <: PrimitiveBlock{2} end

Base.@kwdef struct RydInteract <: AbstractTerm
    atoms::Vector
    C::Real = 2π * 862690
end

Base.@kwdef struct SumOfX <: AbstractTerm
    nsites::Int
    Ω = 1

    function SumOfX(nsites, Ω)
        assert_param(nsites, Ω, :Ω)
        new(nsites, Ω)
    end
end

Base.@kwdef struct SumOfXPhase <: AbstractTerm
    nsites::Int
    Ω = 1
    ϕ

    function SumOfXPhase(nsites, Ω, ϕ)
        assert_param(nsites, Ω, :Ω)
        assert_param(nsites, ϕ, :ϕ)
        new(nsites, Ω, ϕ)
    end
end

Base.@kwdef struct SumOfN <: AbstractTerm
    nsites::Int
    Δ = 1

    function SumOfN(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        new(nsites, Δ)
    end
end

Base.@kwdef struct SumOfZ <: AbstractTerm
    nsites::Int
    Δ = 1

    function SumOfZ(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        new(nsites, Δ)
    end
end

YaoAPI.nqudits(::XPhase) = 1
YaoAPI.nqudits(h::RydInteract) = length(h.atoms)
YaoAPI.nqudits(h::SumOfX) = h.nsites
YaoAPI.nqudits(h::SumOfZ) = h.nsites
YaoAPI.nqudits(h::SumOfN) = h.nsites
