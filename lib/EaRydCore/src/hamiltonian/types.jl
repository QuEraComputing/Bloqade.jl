"""
    AbstractTerm

Abstract term for hamiltonian terms.
"""
abstract type AbstractTerm end

const ConstParamType = Union{Number, AbstractVector{<:Number}, NTuple{N, <:Number} where N}
const ConstParamListType = Union{AbstractVector{<:Number}, NTuple{N, <:Number} where N}

"""
    RydInteract{T<:Number, AtomList <: AbstractVector{<:RydAtom}} <: AbstractTerm
    RydInteract(atoms::AbstractVector{<:RydAtom}, C::Number)

Type for Rydberg interactive term.

# Expression

```math
\\sum_{i, j} \\frac{C}{|r_i - r_j|^6} n_i n_j
```

# Parameters

- `atoms`: a list of atom positions, must be type `RydAtom`, default unit is `μm`.
- `C`: the interaction strength, default unit is `MHz⋅μm^6`. default value is `2π * 109.133 * MHz*µm^6`.
"""
struct RydInteract{T <: Number, AtomList <: AbstractVector{<:RydAtom}} <: AbstractTerm
    atoms::AtomList
    C::T

    function RydInteract{T, AtomList}(atoms::AtomList, C::T) where {T <: Number, AtomList}
        new{T, AtomList}(atoms, C)
    end

    function RydInteract(atoms::AtomList, C) where AtomList
        c = default_unit(MHz * µm^6, C)
        new{typeof(c), AtomList}(atoms, c)
    end
end

RydInteract(atoms::AbstractVector) = RydInteract(atoms, 2π * 109.133 * MHz*µm^6)

"""
    XTerm{Omega, Phi} <: AbstractTerm
    XTerm(nsites::Int, Ωs::Omega, ϕs::Phi)

Type for X term.

# Expression

```math
\\sum_{i} Ω_i (e^{iϕ_i} |0⟩_i⟨1| + e^{-iϕ_i}|1⟩_i⟨0|)
```

# Parameters

- `Ω`: rabi-frequency, the default unit for `Ωs` is `MHz`
- `ϕs` phase, has no unit (or `NoUnits` in `Unitful`).
"""
struct XTerm{Omega, Phi} <: AbstractTerm
    nsites::Int
    Ωs::Omega
    ϕs::Phi

    function XTerm{Omega, Phi}(nsites::Int, Ωs, ϕs) where {Omega, Phi}
        assert_nsites(nsites, Ωs, :Ωs)
        assert_nsites(nsites, ϕs, :ϕs)

        assert_has_time_method(Ωs, "Ωs")
        assert_has_time_method(ϕs, "ϕs")
        new{Omega, Phi}(nsites, Ωs, ϕs)
    end
end

function XTerm(nsites::Int, Ωs, ϕs)
    Ωs = default_unit(MHz, Ωs)
    ϕs = default_unit(NoUnits, ϕs)
    XTerm{typeof(Ωs), typeof(ϕs)}(nsites, Ωs, ϕs)
end

"""
    ZTerm{Delta} <: AbstractTerm
    ZTerm(nsites, Δs::Delta)

Type for Z term.

# Expression

```math
\\sum_i Δ_iσ_i^z
```

# Parameters

- `Δs`: the detuning parameter, the default unit is `MHz`.
"""
struct ZTerm{Delta} <: AbstractTerm
    nsites::Int
    Δs::Delta

    function ZTerm{Delta}(nsites::Int, Δs) where {Delta}
        assert_nsites(nsites, Δs, :Δs)
        new{Delta}(nsites, Δs)
    end
end

function ZTerm(nsites::Int, Δs)
    assert_has_time_method(Δs, "Δs")
    Δs = default_unit(MHz, Δs)
    ZTerm{typeof(Δs)}(nsites, Δs)
end

"""
    NTerm{Delta} <: AbstractTerm
    NTerm(nsites, Δs::Delta)

Type for N term

# Expression

```math
\\sum_i Δ_i n_i
```

# Parameters

- `Δs`: the detuning parameter, the default unit is `MHz`.
"""
struct NTerm{Delta} <: AbstractTerm
    nsites::Int
    Δs::Delta

    function NTerm{Delta}(nsites::Int, Δs) where {Delta}
        assert_nsites(nsites, Δs, :Δs)
        new{Delta}(nsites, Δs)
    end
end

function NTerm(nsites::Int, Δs)
    assert_has_time_method(Δs, "Δs")
    Δs = default_unit(MHz, Δs)
    NTerm{typeof(Δs)}(nsites, Δs)
end

struct Negative{Term} <: AbstractTerm
    term::Term
end

const OrNegative{T} = Union{Negative{T}, T}

Base.:(-)(t::AbstractTerm) = Negative(t)
Base.:(-)(t::Negative) = t.term

struct Hamiltonian{Terms <: Tuple} <: AbstractTerm
    terms::Terms

    function Hamiltonian(terms::Terms) where {Terms}
        first_nsites = nsites(first(terms))
        for idx in 2:length(terms)
            first_nsites == nsites(terms[idx]) ||
                throw(ArgumentError(
                    "nsites mismatch, " *
                    "expect $first_nsites, got $(nsites(terms[idx]))"
                ))
        end
        new{Terms}(terms)
    end
end

# try to infer number of sites from the input
# this is only necessary for CUDA
to_tuple(xs) = (xs..., ) # make it type stable
to_tuple(xs::Tuple) = xs

"""
    XTerm(Ωs::AbstractVector, ϕs::AbstractVector)

Create the `XTerm` from given `Ωs` and `ϕs`.
"""
XTerm(Ωs::AbstractVector, ϕs::AbstractVector) = XTerm(length(Ωs), to_tuple(Ωs), to_tuple(ϕs))

"""
    XTerm(Ωs::Number, ϕs::AbstractVector)

Create the `XTerm` from given `Ωs` and `ϕs`.
"""
XTerm(Ωs, ϕs::AbstractVector) = XTerm(length(ϕs), Ωs, to_tuple(ϕs))

function XTerm(nsites::Int, Ωs::AbstractVector)
    @assert nsites == length(Ωs) "number of sites does not match number of rabi frequencies"
    return XTerm(Ωs)
end

"""
    XTerm(Ωs::AbstractVector, ϕs::Number)

Create the `XTerm` from given `Ωs` and `ϕs`.
"""
XTerm(Ωs::AbstractVector, ϕs) = XTerm(length(Ωs), to_tuple(Ωs), ϕs)

# convenient constructor for simple case
"""
    XTerm(n::Int, Ωs::AbstractVector)

Create a simple `XTerm` from given `Ωs`.
"""
XTerm(Ωs::AbstractVector) = XTerm(length(Ωs), to_tuple(Ωs), nothing)

"""
    XTerm(n::Int, Ω::Number)

Create a simple `XTerm` from given number of atoms `n`
and `Ω`.
"""
XTerm(n::Int, Ω) = XTerm(n, Ω, nothing)

"""
    ZTerm(Δs::AbstractVector)

Create a simple `ZTerm` from given `Δs`.
"""
ZTerm(Δs::AbstractVector) = ZTerm(length(Δs), to_tuple(Δs))

"""
    NTerm(Δs::AbstractVector)

Create a simple `NTerm` from given `Δs`.
"""
NTerm(Δs::AbstractVector) = NTerm(length(Δs), to_tuple(Δs))

"""
    nsites(term)

Return the number of sites of given Hamiltonian term.
"""
function nsites end

nsites(t::XTerm) = t.nsites
nsites(t::ZTerm) = t.nsites
nsites(t::NTerm) = t.nsites
nsites(t::Hamiltonian) = nsites(t.terms[1])
nsites(t::Negative) = nsites(t.term)
nsites(t::RydInteract) = length(t.atoms)

function nsites(terms::Vector{<:AbstractTerm})
    term_nsites = nsites(first(terms))
    for i in 2:length(terms)
        term_nsites == nsites(terms[i]) || error("number of sites is not consistent in the list of hamiltonians")
    end
    return term_nsites
end

hilbert_space(n::Int) = 0:((1<<n)-1)
hilbert_space(t::AbstractTerm) = hilbert_space(nsites(t))

Base.eltype(t::XTerm) = eltype(t.Ωs)
Base.eltype(t::ZTerm) = eltype(t.Δs)
Base.eltype(t::NTerm) = eltype(t.Δs)
Base.eltype(t::Negative) = eltype(t.term)
Base.eltype(t::RydInteract) = typeof(t.C)
Base.eltype(t::Hamiltonian) = promote_type(eltype.(t.terms)...)

function Base.isreal(t::XTerm)
    isnothing(t.ϕs) ? true :
    t.ϕs isa Number ? iszero(t.ϕs) : false
end

Base.isreal(t::ZTerm) = true
Base.isreal(t::NTerm) = true
Base.isreal(t::RydInteract) = true
Base.isreal(t::Negative) = isreal(t.term)
Base.isreal(t::Hamiltonian) = all(isreal, t.terms)

Base.iszero(t::XTerm) = iszero(t.Ωs)
Base.iszero(t::Union{NTerm, ZTerm}) = iszero(t.Δs)
Base.iszero(t::RydInteract) = iszero(t.C)
Base.iszero(t::Negative) = iszero(t.term)
Base.iszero(t::Hamiltonian) = all(iszero, t.terms)

Base.getindex(t::Hamiltonian, i::Int) = t.terms[i]