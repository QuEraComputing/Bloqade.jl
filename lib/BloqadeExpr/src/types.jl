# a linear map for low-level hamiltonian representation

"""
    struct Hamiltonian

`Hamiltonian` stores the dynamic prefactors of each term.
The actual hamiltonian is the sum of `f_i(t) * t_i` where
`f_i` and `t_i` are entries of `fs` and `ts`.
"""
struct Hamiltonian{FS<:Tuple,TS<:Tuple}
    fs::FS # prefactor of each term
    ts::TS # const linear map of each term

    function Hamiltonian(fs, ts)
        all(x -> size(x) == size(ts[1]), ts) || throw(ArgumentError("matrix term should have the same size"))
        return new{typeof(fs),typeof(ts)}(fs, ts)
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
struct StepHamiltonian{T,FS,TS}
    t::T # clock
    h::Hamiltonian{FS,TS}
end

Base.size(h::StepHamiltonian, idx::Int) = size(h.h, idx)
Base.size(h::StepHamiltonian) = size(h.h)

function to_matrix(h::StepHamiltonian)
    return sum(zip(h.h.fs, h.h.ts)) do (f, t)
        return f(h.t) * t
    end
end

function LinearAlgebra.opnorm(h::StepHamiltonian, p = 2)
    return opnorm(to_matrix(h), p)
end

(h::Hamiltonian)(t::Real) = StepHamiltonian(t, h)

# NOTE: these high level expression
# will be converted to lower level Hamiltonian object
# before simulation no need to aggressively specialize on them
# the following operators are mainly used to interface
# with Yao

"""
    XPhase{T, D, name} <: PrimitiveBlock{D}

XPhase operator.

```math
e^{ϕ ⋅ im} |0⟩⟨1| + e^{-ϕ ⋅ im} |1⟩⟨0|
```
"""
struct XPhase{T, D, name} <: PrimitiveBlock{D}
    ϕ::T
end
function XPhase(ϕ::T; nlevel = 2, name = :rydberg) where T
    (name === :rydberg && nlevel in (2, 3)) && return XPhase{T, nlevel, name}(ϕ)
    (name === :hyperfine && nlevel == 3) && return XPhase{T, nlevel, name}(ϕ)
    (name === :hyperfine && nlevel == 2) && throw(ArgumentError("There is no hyperfine operator when `nlevel = 2"))
    throw(ArgumentError("`nlevel` should be 2 or 2 and `name` should be one of `:rydberg` and `:hyperfine`"))
end

struct PuPhase{T, D, name} <: PrimitiveBlock{D}
    ϕ::T
end
function PuPhase(ϕ::T; nlevel = 2, name = :rydberg) where T 
    (name === :rydberg && nlevel in (2, 3)) && return PuPhase{T, nlevel, name}(ϕ)
    (name === :hyperfine && nlevel == 3) && return PuPhase{T, nlevel, name}(ϕ)
    (name === :hyperfine && nlevel == 2) && throw(ArgumentError("There is no hyperfine operator when `nlevel = 2"))
    throw(ArgumentError("`nlevel` should be 2 or 2 and `name` should be one of `:rydberg` and `:hyperfine`"))
end

struct PdPhase{T, D, name} <: PrimitiveBlock{D}
    ϕ::T
end
function PdPhase(ϕ::T; nlevel = 2, name = :rydberg) where T 
    (name === :rydberg && nlevel in (2, 3)) && return PdPhase{T, nlevel, name}(ϕ)
    (name === :hyperfine && nlevel == 3) && return PdPhase{T, nlevel, name}(ϕ)
    (name === :hyperfine && nlevel == 2) && throw(ArgumentError("There is no hyperfine operator when `nlevel = 2"))
    throw(ArgumentError("`nlevel` should be 2 or 2 and `name` should be one of `:rydberg` and `:hyperfine`"))
end

struct X3{name} <: YaoBlocks.ConstantGate{1, 3} end
X3(name) = X3{name}()
struct Z3{name} <: YaoBlocks.ConstantGate{1, 3} end
Z3(name) = Z3{name}()
struct N3{name} <: YaoBlocks.ConstantGate{1, 3} end
N3(name) = N3{name}()
struct Pu3{name} <: YaoBlocks.ConstantGate{1, 3} end
Pu3(name) = Pu3{name}()
struct Pd3{name} <: YaoBlocks.ConstantGate{1, 3} end
Pd3(name) = Pd3{name}()


"""
    AbstractTerm{D} <: PrimitiveBlock{D}

Abstract term for local hamiltonian terms on D-level system.
"""
abstract type AbstractTerm{D} <: PrimitiveBlock{D} end

YaoBlocks.unsafe_getindex(::Type{T}, x::AbstractTerm, i::Integer, j::Integer) where {T,N} = YaoBlocks.unsafe_getindex(T, YaoBlocks.Optimise.to_basictypes(x), i, j)
YaoBlocks.unsafe_getcol(::Type{T}, x::AbstractTerm, j::DitStr{2}) where T = YaoBlocks.unsafe_getcol(T, YaoBlocks.Optimise.to_basictypes(x), j)
YaoBlocks.ishermitian(::AbstractTerm) = true

"""
    struct RydInteract{D} <: AbstractTerm{D}
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
Base.@kwdef struct RydInteract{D} <: AbstractTerm{D}
    atoms::Vector
    C::Real = 2π * 862690

    function RydInteract{D}(atoms, C) where D
        return new{D}(atoms, default_unit(MHz * µm^6, C))
    end
end
RydInteract(atoms, C; nlevel = 2) = RydInteract{nlevel}(atoms, C)

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
Base.@kwdef struct SumOfX{D, name} <: AbstractTerm{D}
    nsites::Int
    Ω = 1

    function SumOfX{D, name}(nsites::Int, Ω) where {D, name}
        assert_param(nsites, Ω, :Ω)
        Ω = default_unit(MHz, Ω)
        return new{D, name}(nsites, Ω)
    end
end
function SumOfX(nsites::Int, Ω; nlevel = 2, name = :rydberg) 
    (name === :rydberg && nlevel in (2, 3)) && return SumOfX{nlevel, name}(nsites, Ω)
    (name === :hyperfine && nlevel == 3) && return SumOfX{nlevel, name}(nsites, Ω)
    (name === :hyperfine && nlevel == 2) && throw(ArgumentError("There is no hyperfine operator when `nlevel = 2"))
    throw(ArgumentError("`nlevel` should be 2 or 2 and `name` should be one of `:rydberg` and `:hyperfine`"))
end

SumOfX(n::Int; nlevel = 2, name = :rydberg) = SumOfX(n, 1; nlevel = nlevel, name = name)

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
Base.@kwdef struct SumOfXPhase{D, name} <: AbstractTerm{D}
    nsites::Int
    Ω = 1
    ϕ::Any

    function SumOfXPhase{D, name}(nsites, Ω, ϕ) where {D, name}
        assert_param(nsites, Ω, :Ω)
        assert_param(nsites, ϕ, :ϕ)

        Ω = default_unit(MHz, Ω)
        ϕ = default_unit(NoUnits, ϕ)
        return new{D, name}(nsites, Ω, ϕ)
    end
end
function SumOfXPhase(nsites::Int, Ω, ϕ; nlevel = 2, name = :rydberg) 
    (name === :rydberg && nlevel in (2, 3)) && return SumOfXPhase{nlevel, name}(nsites, Ω, ϕ)
    (name === :hyperfine && nlevel == 3) && return SumOfXPhase{nlevel, name}(nsites, Ω, ϕ)
    (name === :hyperfine && nlevel == 2) && throw(ArgumentError("There is no hyperfine operator when `nlevel = 2"))
    throw(ArgumentError("`nlevel` should be 2 or 2 and `name` should be one of `:rydberg` and `:hyperfine`"))
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
Base.@kwdef struct SumOfN{D, name} <: AbstractTerm{D}
    nsites::Int
    Δ = 1

    function SumOfN{D, name}(nsites, Δ) where {D, name}
        assert_param(nsites, Δ, :Δ)
        return new{D, name}(nsites, default_unit(MHz, Δ))
    end
end

function SumOfN(nsites::Int, Δ; nlevel = 2, name = :rydberg) 
    (name === :rydberg && nlevel in (2, 3)) && return SumOfN{nlevel, name}(nsites, Δ)
    (name === :hyperfine && nlevel == 3) && return SumOfN{nlevel, name}(nsites, Δ)
    (name === :hyperfine && nlevel == 2) && throw(ArgumentError("There is no hyperfine operator when `nlevel = 2"))
    throw(ArgumentError("`nlevel` should be 2 or 2 and `name` should be one of `:rydberg` and `:hyperfine`"))
end
SumOfN(n::Int; nlevel = 2, name = :rydberg) = SumOfN(n, 1; nlevel = nlevel, name = name)

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
Base.@kwdef struct SumOfZ{D, name} <: AbstractTerm{D}
    nsites::Int
    Δ = 1

    function SumOfZ{D, name}(nsites, Δ) where {D, name}
        assert_param(nsites, Δ, :Δ)
        return new{D, name}(nsites, default_unit(MHz, Δ))
    end
end

function SumOfZ(nsites::Int, Δ; nlevel = 2, name = :rydberg) 
    (name === :rydberg && nlevel in (2, 3)) && return SumOfZ{nlevel, name}(nsites, Δ)
    (name === :hyperfine && nlevel == 3) && return SumOfZ{nlevel, name}(nsites, Δ)
    (name === :hyperfine && nlevel == 2) && throw(ArgumentError("There is no hyperfine operator when `nlevel = 2"))
    throw(ArgumentError("`nlevel` should be 2 or 2 and `name` should be one of `:rydberg` and `:hyperfine`"))
end
SumOfZ(n::Int; nlevel = 2, name = :rydberg) = SumOfZ(n, 1; nlevel = nlevel, name = name)

struct DivByTwo{F}
    f::F
    function DivByTwo(f)
        assert_has_time_method(f,:f)
        return new{typeof(f)}(f)
    end
end

@inline (Func::DivByTwo)(t::Real) = Func.f(t)/2


const RabiTypes = Union{Nothing,SumOfX,SumOfXPhase}
const DetuningTypes = Union{Nothing,SumOfN}

struct RydbergHamiltonian{D} <: AbstractTerm{D}
    rydberg_term::RydInteract{D}
    rabi_term::RabiTypes
    detuning_term::DetuningTypes
end

function add_terms(h::RydbergHamiltonian)
    terms = h.rydberg_term

    if typeof(h.rabi_term) <: Union{SumOfX,SumOfXPhase}
        terms += h.rabi_term
    end
    
    if h.detuning_term isa SumOfN
        terms -= h.detuning_term
    end
    
    return YaoBlocks.Optimise.simplify(terms)

end

function YaoBlocks.unsafe_getindex(::Type{T}, h::RydbergHamiltonian, i::Integer, j::Integer) where {T,N}
    return YaoBlocks.unsafe_getindex(T, YaoBlocks.Optimise.to_basictypes(h), i, j)
end

function YaoBlocks.unsafe_getcol(::Type{T}, h::RydbergHamiltonian, j::DitStr{2}) where T
    YaoBlocks.unsafe_getcol(T, YaoBlocks.Optimise.to_basictypes(h), j)
end

YaoAPI.nqudits(::XPhase) = 1
YaoAPI.nqudits(::PdPhase) = 1
YaoAPI.nqudits(::PuPhase) = 1
YaoAPI.nqudits(h::RydInteract) = length(h.atoms)
YaoAPI.nqudits(h::SumOfX) = h.nsites
YaoAPI.nqudits(h::SumOfXPhase) = h.nsites
YaoAPI.nqudits(h::SumOfZ) = h.nsites
YaoAPI.nqudits(h::SumOfN) = h.nsites
@inline YaoAPI.nqudits(h::RydbergHamiltonian) = nqudits(h.rydberg_term)


function Base.:(==)(lhs::RydInteract, rhs::RydInteract)
    return lhs.C == rhs.C && lhs.atoms == rhs.atoms
end

function Base.:(==)(lhs::SumOfX, rhs::SumOfX)
    return lhs.nsites == rhs.nsites && lhs.Ω == rhs.Ω
end

function Base.:(==)(lhs::SumOfZ, rhs::SumOfZ)
    return lhs.nsites == rhs.nsites && lhs.Δ == rhs.Δ
end

function Base.:(==)(lhs::SumOfN, rhs::SumOfN)
    return lhs.nsites == rhs.nsites && lhs.Δ == rhs.Δ
end

function Base.:(==)(lhs::SumOfXPhase, rhs::SumOfXPhase)
    return lhs.nsites == rhs.nsites && lhs.Ω == rhs.Ω && lhs.ϕ == rhs.ϕ
end

function Base.:(==)(lhs::RydbergHamiltonian, rhs::RydbergHamiltonian)
    return lhs.rydberg_term == rhs.rydberg_term && lhs.rabi_term == rhs.rabi_term && lhs.detuning_term == rhs.detuning_term
end

Base.isreal(::RydInteract) = true
Base.isreal(::SumOfN) = true
Base.isreal(::SumOfX) = true
Base.isreal(::SumOfZ) = true
Base.isreal(::SumOfXPhase) = false
Base.isreal(h::Add) = all(isreal, subblocks(h))
Base.isreal(h::Scale) = isreal(factor(h)) && isreal(content(h))
Base.isreal(h::RydbergHamiltonian) = !(h.rabi_term isa SumOfXPhase)



storage_size(x) = sizeof(x)
function storage_size(h::Hamiltonian)
    return sum(storage_size, h.ts)
end
function storage_size(H::SparseMatrixCSC)
    return sizeof(H.colptr) + sizeof(H.rowval) + sizeof(H.nzval)
end
