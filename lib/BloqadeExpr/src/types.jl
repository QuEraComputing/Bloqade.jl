# a linear map for low-level hamiltonian representation

"""
    struct Hamiltonian

`Hamiltonian` stores the dynamic prefactors of each term.
The actual hamiltonian is the sum of `f_i(t) * t_i` where
`f_i` and `t_i` are entries of `fs` and `ts`.
"""
struct Hamiltonian{FS <: Tuple, TS <: Tuple}
    fs::FS # prefactor of each term
    ts::TS # const linear map of each term

    function Hamiltonian(fs, ts)
        all(x->size(x)==size(ts[1]), ts) || throw(ArgumentError("matrix term should have the same size"))
        new{typeof(fs), typeof(ts)}(fs, ts)
    end
end

Base.size(h::Hamiltonian) = size(h.ts[1])
Base.size(h::Hamiltonian, idx::Int) = size(h.ts[1], idx)

Adapt.@adapt_structure Hamiltonian

"""
    struct StepHamiltonian

A low-level linear-map object that encodes time-dependent
hamiltonian at time step `t`. This object supports the
linear map interface `mul!(Y, H, X)`.
"""
struct StepHamiltonian{T, FS, TS}
    t::T # clock
    h::Hamiltonian{FS, TS}
end

Base.size(h::StepHamiltonian, idx::Int) = size(h.h, idx)
Base.size(h::StepHamiltonian) = size(h.h)

function to_matrix(h::StepHamiltonian)
    return sum(zip(h.h.fs, h.h.ts)) do (f, t)
        f(h.t) * t
    end
end

function LinearAlgebra.opnorm(h::StepHamiltonian, p=2)
    return opnorm(to_matrix(h), p)
end

(h::Hamiltonian)(t::Real) = StepHamiltonian(t, h)

# NOTE: these high level expression
# will be converted to lower level Hamiltonian object
# before simulation no need to aggressively specialize on them
# the following operators are mainly used to interface
# with Yao

"""
    XPhase{T} <: PrimitiveBlock{2}

XPhase operator.

```math
e^{ϕ ⋅ im} |0⟩⟨1| + e^{-ϕ ⋅ im} |1⟩⟨0|
```
"""
struct XPhase{T} <: PrimitiveBlock{2}
    ϕ::T
end

struct PuPhase{T} <: PrimitiveBlock{2}
    ϕ::T
end

struct PdPhase{T} <: PrimitiveBlock{2}
    ϕ::T
end



"""
    AbstractTerm <: PrimitiveBlock{2}

Abstract term for local hamiltonian terms.
"""
abstract type AbstractTerm <: PrimitiveBlock{2} end

"""
    struct RydInteract <: AbstractTerm
    RydInteract(;atoms, C=2π * 862690MHz⋅μm^6)

Type for Rydberg interactive term.

# Expression

```math
\\sum_{i, j} \\frac{C}{|r_i - r_j|^6} n_i n_j
```

# Keyword Arguments

- `atoms`: a list of atom positions, must be type `RydAtom`, default unit is `μm`.
- `C`: the interaction strength, default unit is `MHz⋅μm^6`. default value is `2π * 862690 * MHz*µm^6`.
"""
Base.@kwdef struct RydInteract <: AbstractTerm
    atoms::Vector
    C::Real = 2π * 862690

    function RydInteract(atoms, C)
        new(atoms, default_unit(MHz * µm^6, C))
    end
end

"""
    struct SumOfX <: AbstractTerm
    SumOfX(nsites, Ω)

Term for sum of X operators.

The following two expressions are equivalent

```jldoctest; setup=:(using BloqadeExpr, YaoBlocks)
julia> SumOfX(nsites=5)
∑ σ^x_i

julia> sum([X for _ in 1:5])
nqudits: 1
+
├─ X
├─ X
├─ X
├─ X
└─ X
```

# Expression

```math
\\sum_i Ω σ^x_i
```
"""
Base.@kwdef struct SumOfX <: AbstractTerm
    nsites::Int
    Ω = 1

    function SumOfX(nsites, Ω)
        assert_param(nsites, Ω, :Ω)
        Ω = default_unit(MHz, Ω)
        new(nsites, Ω)
    end
end

SumOfX(n::Int) = SumOfX(n, 1)

"""
    struct SumOfXPhase <: AbstractTerm
    SumOfXPhase(;nsites, Ω=1, ϕ)

Sum of `XPhase` operators.

The following two expressions are equivalent

```jldoctest; setup=:(using BloqadeExpr)
julia> SumOfXPhase(nsites=5, ϕ=0.1)
1.0 ⋅ ∑ e^{0.1 ⋅ im} |0⟩⟨1| + e^{-0.1 ⋅ im} |1⟩⟨0|

julia> sum([XPhase(0.1) for _ in 1:5])
nqudits: 1
+
├─ XPhase(0.1)
├─ XPhase(0.1)
├─ XPhase(0.1)
├─ XPhase(0.1)
└─ XPhase(0.1)
```

But may provide extra speed up.

# Expression

```math
\\sum_i Ω ⋅ (e^{ϕ ⋅ im} |0⟩⟨1| + e^{-ϕ ⋅ im} |1⟩⟨0|)
```
"""
Base.@kwdef struct SumOfXPhase <: AbstractTerm
    nsites::Int
    Ω = 1
    ϕ

    function SumOfXPhase(nsites, Ω, ϕ)
        assert_param(nsites, Ω, :Ω)
        assert_param(nsites, ϕ, :ϕ)

        Ω = default_unit(MHz, Ω)
        ϕ = default_unit(NoUnits, ϕ)
        new(nsites, Ω, ϕ)
    end
end

"""
    struct SumOfN <: AbstractTerm
    SumOfN(;nsites[, Δ=1])

Sum of N operators. 

The following two expression are equivalent

```jldoctest; setup=:(using BloqadeExpr)
julia> SumOfN(nsites=5)
∑ n_i

julia> sum([Op.n for _ in 1:5])
nqudits: 1
+
├─ P1
├─ P1
├─ P1
├─ P1
└─ P1
```

But may provide extra speed up.

# Expression

```math
\\sum_i Δ ⋅ n_i
```
"""
Base.@kwdef struct SumOfN <: AbstractTerm
    nsites::Int
    Δ = 1

    function SumOfN(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        new(nsites, default_unit(MHz, Δ))
    end
end

SumOfN(n::Int) = SumOfN(n, 1)

"""
    struct SumOfZ <: AbstractTerm
    SumOfZ(;nsites, Δ=1)

Sum of Pauli Z operators.

The following two expression are equivalent

```jldoctest; setup=:(using BloqadeExpr, YaoBlocks)
julia> SumOfZ(nsites=5)
∑ σ^z_i

julia> sum([Z for _ in 1:5])
nqudits: 1
+
├─ Z
├─ Z
├─ Z
├─ Z
└─ Z
```

# Expression

```math
\\sum_i Δ ⋅ σ^z_i
```
"""
Base.@kwdef struct SumOfZ <: AbstractTerm
    nsites::Int
    Δ = 1

    function SumOfZ(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        new(nsites, default_unit(MHz, Δ))
    end
end

SumOfZ(n::Int) = SumOfZ(n, 1)

YaoAPI.nqudits(::XPhase) = 1
YaoAPI.nqudits(::PdPhase) = 1
YaoAPI.nqudits(::PuPhase) = 1
YaoAPI.nqudits(h::RydInteract) = length(h.atoms)
YaoAPI.nqudits(h::SumOfX) = h.nsites
YaoAPI.nqudits(h::SumOfXPhase) = h.nsites
YaoAPI.nqudits(h::SumOfZ) = h.nsites
YaoAPI.nqudits(h::SumOfN) = h.nsites

function Base.:(==)(lhs::RydInteract, rhs::RydInteract)
    lhs.C == rhs.C && lhs.atoms == rhs.atoms
end

function Base.:(==)(lhs::SumOfX, rhs::SumOfX)
    lhs.nsites == rhs.nsites && lhs.Ω == rhs.Ω
end

function Base.:(==)(lhs::SumOfZ, rhs::SumOfZ)
    lhs.nsites == rhs.nsites && lhs.Δ == rhs.Δ
end

function Base.:(==)(lhs::SumOfN, rhs::SumOfN)
    lhs.nsites == rhs.nsites && lhs.Δ == rhs.Δ
end

function Base.:(==)(lhs::SumOfXPhase, rhs::SumOfXPhase)
    lhs.nsites == rhs.nsites && lhs.Ω == rhs.Ω && lhs.ϕ == rhs.ϕ
end

Base.isreal(::RydInteract) = true
Base.isreal(::SumOfN) = true
Base.isreal(::SumOfX) = true
Base.isreal(::SumOfZ) = true
Base.isreal(::SumOfXPhase) = false
Base.isreal(h::Add) = all(isreal, subblocks(h))
Base.isreal(h::Scale) = isreal(factor(h)) && isreal(content(h))

storage_size(x) = sizeof(x)
function storage_size(h::Hamiltonian)
    sum(storage_size, h.ts)
end
function storage_size(H::SparseMatrixCSC)
    return sizeof(H.colptr) + sizeof(H.rowval) + sizeof(H.nzval)
end
