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
    XPhase{T} <: PrimitiveBlock{2}

XPhase operator for 2-level Rydberg system.

```math
e^{ϕ ⋅ i} |0⟩⟨1| + e^{-ϕ ⋅ i} |1⟩⟨0|
```
"""
struct XPhase{T} <: PrimitiveBlock{2}
    ϕ::T
end

"""
    XPhase_01{T} <: PrimitiveBlock{3}

XPhase operator act on |0⟩ and |1⟩ for 3-level Rydberg system.

```math
e^{ϕ ⋅ i} |0⟩⟨1| + e^{-ϕ ⋅ i} |1⟩⟨0| = 
\\begin{pmatrix}
0 & e^{ϕ ⋅ im} & 0 \\\\
e^{-ϕ ⋅ im} & 0 & 0 \\\\
0 & 0 & 0
\\end{pmatrix}
```
"""
struct XPhase_01{T} <: PrimitiveBlock{3}
    ϕ::T
end

"""
    XPhase_1r{T} <: PrimitiveBlock{3}

XPhase operator act on |1⟩ and |r⟩ for 3-level Rydberg system.

```math
e^{ϕ ⋅ i} |1⟩⟨r| + e^{-ϕ ⋅ i} |r⟩⟨1| = 
\\begin{pmatrix}
0 & 0 & 0 \\\\
0 & 0 & e^{ϕ ⋅ im} \\\\
0 & e^{-ϕ ⋅ im} & 0
\\end{pmatrix}
```
"""
struct XPhase_1r{T} <: PrimitiveBlock{3}
    ϕ::T
end

struct PuPhase{T} <: PrimitiveBlock{2}
    ϕ::T
end
struct PuPhase_01{T} <: PrimitiveBlock{3}
    ϕ::T
end
struct PuPhase_1r{T} <: PrimitiveBlock{3}
    ϕ::T
end

struct PdPhase{T} <: PrimitiveBlock{2}
    ϕ::T
end
struct PdPhase_01{T} <: PrimitiveBlock{3}
    ϕ::T
end
struct PdPhase_1r{T} <: PrimitiveBlock{3}
    ϕ::T
end

struct OpX_01 <: YaoBlocks.ConstantGate{1, 3} end
struct OpX_1r <: YaoBlocks.ConstantGate{1, 3} end
struct OpZ_01 <: YaoBlocks.ConstantGate{1, 3} end
struct OpZ_1r <: YaoBlocks.ConstantGate{1, 3} end
struct OpN_1 <: YaoBlocks.ConstantGate{1, 3} end
struct OpN_r <: YaoBlocks.ConstantGate{1, 3} end
struct OpPu_01 <: YaoBlocks.ConstantGate{1, 3} end
struct OpPu_1r <: YaoBlocks.ConstantGate{1, 3} end
struct OpPd_01 <: YaoBlocks.ConstantGate{1, 3} end
struct OpPd_1r <: YaoBlocks.ConstantGate{1, 3} end

const X_01 = OpX_01()
const X_1r = OpX_1r()
const Z_01 = OpZ_01()
const Z_1r = OpZ_1r()
const N_1 = OpN_1()
const N_r = OpN_r()
const Pu_01 = OpPu_01()
const Pu_1r = OpPu_1r()
const Pd_01 = OpPd_01()
const Pd_1r = OpPd_1r()


# ==================== Docstring for 3-level constant gates ====================

"""
    X_01
    OpX_01 <: YaoBlocks.ConstantGate{1, 3}

Pauli X operator act on |0⟩ and |1⟩ for 3-level Rydberg system.

Matrix expression:

```math
\\sigma^{x,\\mathrm{hf}} = 
\\begin{pmatrix}
0 & 1 & 0 \\\\
1 & 0 & 0 \\\\
0 & 0 & 0
\\end{pmatrix}
```
"""
X_01, OpX_01

"""
    X_1r
    OpX_1r <: YaoBlocks.ConstantGate{1, 3}

Pauli X operator act on |1⟩ and |r⟩ for 3-level Rydberg system.

Matrix expression:

```math
\\sigma^{x,\\mathrm{r}} = 
\\begin{pmatrix}
0 & 0 & 0 \\\\
0 & 0 & 1 \\\\
0 & 1 & 0
\\end{pmatrix}
```
"""
X_1r, OpX_1r

"""
    Z_01
    OpZ_01 <: YaoBlocks.ConstantGate{1, 3}

Pauli Z operator act on |0⟩ and |1⟩ for 3-level Rydberg system.

Matrix expression:

```math
\\sigma^{z,\\mathrm{hf}} = 
\\begin{pmatrix}
1 & 0 & 0 \\\\
0 & -1 & 0 \\\\
0 & 0 & 0
\\end{pmatrix}
```
"""
Z_01, OpZ_01

"""
    Z_1r
    OpZ_1r <: YaoBlocks.ConstantGate{1, 3}

Pauli Z operator act on |1⟩ and |r⟩ for 3-level Rydberg system.

Matrix expression:

```math
\\sigma^{z,\\mathrm{r}} = 
\\begin{pmatrix}
0 & 0 & 0 \\\\
0 & 1 & 0 \\\\
0 & 0 & -1
\\end{pmatrix}
```
"""
Z_1r, OpZ_1r

"""
    N_1
    OpN_1 <: YaoBlocks.ConstantGate{1, 3}

Projection operator onto |1⟩ for 3-level Rydberg system.

Matrix expression:

```math
n^1 = |1⟩⟨1| = 
\\begin{pmatrix}
0 & 0 & 0 \\\\
0 & 1 & 0 \\\\
0 & 0 & 0
\\end{pmatrix}
```
"""
N_1, OpN_1

"""
    N_r
    OpN_r <: YaoBlocks.ConstantGate{1, 3}

Projection operator onto |r⟩ for 3-level Rydberg system.

Matrix expression:
    
```math
n^{\\mathrm{r}} = |r⟩⟨r| = 
\\begin{pmatrix}
0 & 0 & 0 \\\\
0 & 0 & 0 \\\\
0 & 0 & 1
\\end{pmatrix}
```
"""
N_r, OpN_r

"""
    Pu_01
    OpPu_01 <: YaoBlocks.ConstantGate{1, 3}

Matrix expression:

```math
\\mathrm{Pu}^{\\mathrm{hf}} = 
\\begin{pmatrix}
0 & 1 & 0 \\\\
0 & 0 & 0 \\\\
0 & 0 & 0
\\end{pmatrix}
```
"""
Pu_01, OpPu_01


"""
    Pu_1r
    OpPu_1r <: YaoBlocks.ConstantGate{1, 3}

Matrix expression:

```math
\\mathrm{Pu}^{\\mathrm{r}} = 
\\begin{pmatrix}
0 & 0 & 0 \\\\
0 & 0 & 1 \\\\
0 & 0 & 0
\\end{pmatrix}
```
"""
Pu_1r, OpPu_1r

"""
    Pd_01
    OpPd_01 <: YaoBlocks.ConstantGate{1, 3}

Matrix expression:

```math
\\mathrm{Pd}^{\\mathrm{hf}} = 
\\begin{pmatrix}
0 & 0 & 0 \\\\
1 & 0 & 0 \\\\
0 & 0 & 0
\\end{pmatrix}
```
"""
Pd_01, OpPd_01

"""
    Pd_1r
    OpPd_1r <: YaoBlocks.ConstantGate{1, 3}

Matrix expression:

```math
\\mathrm{Pd}^{\\mathrm{r}} = 
\\begin{pmatrix}
0 & 0 & 0 \\\\
0 & 0 & 0 \\\\
0 & 1 & 0
\\end{pmatrix}
```
"""
Pd_1r, OpPd_1r

# ==================== Docstring for 3-level constant gates (end) ====================


"""
    AbstractTerm{D} <: PrimitiveBlock{D}

Abstract term for local hamiltonian terms on D-level system.
"""
abstract type AbstractTerm{D} <: PrimitiveBlock{D} end

YaoBlocks.unsafe_getindex(::Type{T}, x::AbstractTerm, i::Integer, j::Integer) where T = YaoBlocks.unsafe_getindex(T, YaoBlocks.Optimise.to_basictypes(x), i, j)
YaoBlocks.unsafe_getcol(::Type{T}, x::AbstractTerm, j::DitStr{2}) where T = YaoBlocks.unsafe_getcol(T, YaoBlocks.Optimise.to_basictypes(x), j)
YaoBlocks.ishermitian(::AbstractTerm) = true

"""
    struct RydInteract{D} <: AbstractTerm{D}
    RydInteract(;atoms, C=2π * 862690MHz⋅μm^6)

Type for Rydberg interactive term.

# Expression

```math
\\sum_{i, j} \\frac{C}{|x_i - x_j|^6} n_i n_j
```

# Keyword Arguments

- `atoms`: a list of atom positions, must be type `RydAtom`, default unit is `μm`.
- `C`: the interaction strength, default unit is `MHz⋅μm^6`. default value is `2π * 862690 * MHz*µm^6`.
"""
Base.@kwdef struct RydInteract{D} <: AbstractTerm{D}
    atoms::Union{Vector,BoundedLattice}
    C::Real = 2π * 862690

    function RydInteract{D}(atoms, C) where D
        return new{D}(atoms, default_unit(MHz * µm^6, C))
    end
end
RydInteract(atoms, C; nlevel = 2) = RydInteract{nlevel}(atoms, C)

"""
    struct SumOfX <: AbstractTerm{2}
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
Base.@kwdef struct SumOfX <: AbstractTerm{2}
    nsites::Int
    Ω = 1

    function SumOfX(nsites::Int, Ω)
        assert_param(nsites, Ω, :Ω)
        Ω = default_unit(MHz, Ω)
        return new(nsites, Ω)
    end
end
SumOfX(n::Int) = SumOfX(n, 1)

"""
    struct SumOfX_01 <: AbstractTerm{3}
    SumOfX_01(nsites, Ω)

Term for sum of `X_01` operators.

# Expression

```math
\\sum_i Ω σ^{x,\\mathrm{hf}}_i
```
"""
Base.@kwdef struct SumOfX_01 <: AbstractTerm{3}
    nsites::Int
    Ω = 1

    function SumOfX_01(nsites::Int, Ω)
        assert_param(nsites, Ω, :Ω)
        Ω = default_unit(MHz, Ω)
        return new(nsites, Ω)
    end
end
SumOfX_01(n::Int) = SumOfX_01(n, 1)

"""
    struct SumOfX_1r <: AbstractTerm{3}
    SumOfX_1r(nsites, Ω)

Term for sum of `X_1r` operators.

# Expression

```math
\\sum_i Ω σ^{x,\\mathrm{r}}_i
```
"""
Base.@kwdef struct SumOfX_1r <: AbstractTerm{3}
    nsites::Int
    Ω = 1

    function SumOfX_1r(nsites::Int, Ω)
        assert_param(nsites, Ω, :Ω)
        Ω = default_unit(MHz, Ω)
        return new(nsites, Ω)
    end
end
SumOfX_1r(n::Int) = SumOfX_1r(n, 1)


"""
    struct SumOfXPhase <: AbstractTerm{2}
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
\\sum_i Ω ⋅ (e^{ϕ ⋅ i} |0⟩⟨1| + e^{-ϕ ⋅ i} |1⟩⟨0|)
```
"""
Base.@kwdef struct SumOfXPhase <: AbstractTerm{2}
    nsites::Int
    Ω = 1
    ϕ::Any

    function SumOfXPhase(nsites, Ω, ϕ)
        assert_param(nsites, Ω, :Ω)
        assert_param(nsites, ϕ, :ϕ)

        Ω = default_unit(MHz, Ω)
        ϕ = default_unit(NoUnits, ϕ)
        return new(nsites, Ω, ϕ)
    end
end

"""
    struct SumOfXPhase_01 <: AbstractTerm{3}
    SumOfXPhase_01(nsites, Ω, ϕ)

Term for sum of `XPhase_01` operators.

# Expression

```math
\\sum_i Ω ⋅ (e^{ϕ ⋅ i} |0⟩⟨1| + e^{-ϕ ⋅ i} |1⟩⟨0|)
```
"""
Base.@kwdef struct SumOfXPhase_01 <: AbstractTerm{3}
    nsites::Int
    Ω = 1
    ϕ::Any

    function SumOfXPhase_01(nsites, Ω, ϕ)
        assert_param(nsites, Ω, :Ω)
        assert_param(nsites, ϕ, :ϕ)

        Ω = default_unit(MHz, Ω)
        ϕ = default_unit(NoUnits, ϕ)
        return new(nsites, Ω, ϕ)
    end
end

"""
    struct SumOfXPhase_1r <: AbstractTerm{3}
    SumOfXPhase_1r(nsites, Ω, ϕ)

Term for sum of `XPhase_1r` operators.

# Expression

```math
\\sum_i Ω ⋅ (e^{ϕ ⋅ i} |1⟩⟨r| + e^{-ϕ ⋅ i} |r⟩⟨1|)
```
"""
Base.@kwdef struct SumOfXPhase_1r <: AbstractTerm{3}
    nsites::Int
    Ω = 1
    ϕ::Any

    function SumOfXPhase_1r(nsites, Ω, ϕ)
        assert_param(nsites, Ω, :Ω)
        assert_param(nsites, ϕ, :ϕ)

        Ω = default_unit(MHz, Ω)
        ϕ = default_unit(NoUnits, ϕ)
        return new(nsites, Ω, ϕ)
    end
end


"""
    struct SumOfN <: AbstractTerm{2}
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
Base.@kwdef struct SumOfN <: AbstractTerm{2}
    nsites::Int
    Δ = 1

    function SumOfN(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        return new(nsites, default_unit(MHz, Δ))
    end
end
SumOfN(n::Int) = SumOfN(n, 1)

"""
    struct SumOfN_r <: AbstractTerm{3}
    SumOfN_1(;nsites[, Δ=1])

Sum of N_r operators. 

# Expression

```math
\\sum_i Δ ⋅ n^r_i
```
"""
Base.@kwdef struct SumOfN_r <: AbstractTerm{3}
    nsites::Int
    Δ = 1

    function SumOfN_r(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        return new(nsites, default_unit(MHz, Δ))
    end
end
SumOfN_r(n::Int) = SumOfN_r(n, 1)

"""
    struct SumOfN_1 <: AbstractTerm{3}
    SumOfN_1(;nsites[, Δ=1])

Sum of N_1 operators. 

# Expression

```math
\\sum_i Δ ⋅ n^r_i
```
"""
Base.@kwdef struct SumOfN_1 <: AbstractTerm{3}
    nsites::Int
    Δ = 1

    function SumOfN_1(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        return new(nsites, default_unit(MHz, Δ))
    end
end
SumOfN_1(n::Int) = SumOfN_1(n, 1)

"""
    struct SumOfZ <: AbstractTerm{2}
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
Base.@kwdef struct SumOfZ <: AbstractTerm{2}
    nsites::Int
    Δ = 1

    function SumOfZ(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        return new(nsites, default_unit(MHz, Δ))
    end
end
SumOfZ(n::Int) = SumOfZ(n, 1)

"""
    struct SumOfZ_01 <: AbstractTerm{2}
    SumOfZ_01(;nsites, Δ=1)

Sum of Pauli Z_01 operators.

# Expression

```math
\\sum_i Δ ⋅ σ^{z,\\mathrm{hf}}_i
```
"""
Base.@kwdef struct SumOfZ_01 <: AbstractTerm{3}
    nsites::Int
    Δ = 1

    function SumOfZ_01(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        return new(nsites, default_unit(MHz, Δ))
    end
end
SumOfZ_01(n::Int) = SumOfZ_01(n, 1)

"""
    struct SumOfZ_1r <: AbstractTerm{2}
    SumOfZ_1r(;nsites, Δ=1)

Sum of Pauli Z_1r operators.

# Expression

```math
\\sum_i Δ ⋅ σ^{z,\\mathrm{r}}_i
```
"""
Base.@kwdef struct SumOfZ_1r <: AbstractTerm{3}
    nsites::Int
    Δ = 1

    function SumOfZ_1r(nsites, Δ)
        assert_param(nsites, Δ, :Δ)
        return new(nsites, default_unit(MHz, Δ))
    end
end
SumOfZ_1r(n::Int) = SumOfZ_1r(n, 1)


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

const HyperfineRabiTypes = Union{Nothing,SumOfX_01,SumOfXPhase_01}
const HyperfineDetuningTypes = Union{Nothing,SumOfN_1}
const RydbergRabiTypes = Union{Nothing,SumOfX_1r,SumOfXPhase_1r}
const RydbergDetuningTypes = Union{Nothing,SumOfN_r}

const SumOfXPhaseTypes = Union{SumOfXPhase, SumOfXPhase_01, SumOfXPhase_1r}
const SumOfZAndNTypes = Union{SumOfZ, SumOfZ_01, SumOfZ_1r, SumOfN, SumOfN_1, SumOfN_r}
const SumOfXTypes = Union{SumOfX, SumOfX_01, SumOfX_1r}
const ThreeLevelRydbergConstGates = Union{OpX_01, OpX_1r, OpZ_01, OpZ_1r, OpN_1, OpN_r, OpPu_01, OpPu_1r, OpPd_01, OpPd_1r}


struct RydbergHamiltonian <: AbstractTerm{2}
    rydberg_term::RydInteract{2}
    rabi_term::RabiTypes
    detuning_term::DetuningTypes
end

struct RydbergHamiltonian3 <: AbstractTerm{3}
    rydberg_term::RydInteract{3}
    rabi_term_hf::HyperfineRabiTypes
    detuning_term_hf::HyperfineDetuningTypes
    rabi_term_r::RydbergRabiTypes
    detuning_term_r::RydbergDetuningTypes
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

function add_terms(h::RydbergHamiltonian3)
    terms = h.rydberg_term

    if typeof(h.rabi_term_hf) <: Union{SumOfX_01,SumOfXPhase_01}
        terms += h.rabi_term_hf
    end
    if typeof(h.rabi_term_r) <: Union{SumOfX_1r,SumOfXPhase_1r}
        terms += h.rabi_term_r
    end
    if h.detuning_term_hf isa SumOfN_1
        terms -= h.detuning_term_hf
    end
    if h.detuning_term_r isa SumOfN_r
        terms -= h.detuning_term_r
    end
    
    return YaoBlocks.Optimise.simplify(terms)
end

function YaoBlocks.unsafe_getindex(::Type{T}, h::RydbergHamiltonian, i::Integer, j::Integer) where T
    return YaoBlocks.unsafe_getindex(T, YaoBlocks.Optimise.to_basictypes(h), i, j)
end

function YaoBlocks.unsafe_getcol(::Type{T}, h::RydbergHamiltonian, j::DitStr{2}) where T
    YaoBlocks.unsafe_getcol(T, YaoBlocks.Optimise.to_basictypes(h), j)
end


YaoAPI.nqudits(::XPhase) = 1
YaoAPI.nqudits(::XPhase_01) = 1
YaoAPI.nqudits(::XPhase_1r) = 1
YaoAPI.nqudits(::PdPhase) = 1
YaoAPI.nqudits(::PdPhase_01) = 1
YaoAPI.nqudits(::PdPhase_1r) = 1
YaoAPI.nqudits(::PuPhase) = 1
YaoAPI.nqudits(::PuPhase_01) = 1
YaoAPI.nqudits(::PuPhase_1r) = 1
YaoAPI.nqudits(h::RydInteract) = length(h.atoms)
YaoAPI.nqudits(h::SumOfXTypes) = h.nsites
YaoAPI.nqudits(h::SumOfXPhaseTypes) = h.nsites
YaoAPI.nqudits(h::SumOfZAndNTypes) = h.nsites
YaoAPI.nqudits(::ThreeLevelRydbergConstGates) = 1
@inline YaoAPI.nqudits(h::Union{RydbergHamiltonian, RydbergHamiltonian3}) = nqudits(h.rydberg_term)

# checks of objects have the same base object f.
function Base.:(==)(lhs::DivByTwo{F1},rhs::DivByTwo{F2}) where {F1,F2} 
    return F1 == F2 && lhs.f == rhs.f
end

function Base.:(==)(lhs::RydInteract{D1}, rhs::RydInteract{D2}) where {D1, D2}
    return lhs.C == rhs.C && lhs.atoms == rhs.atoms && D1 == D2
end

function Base.:(==)(lhs::SumOfXTypes, rhs::SumOfXTypes)
    return lhs.nsites == rhs.nsites && lhs.Ω == rhs.Ω && typeof(lhs) == typeof(rhs)
end

function Base.:(==)(lhs::SumOfZAndNTypes, rhs::SumOfZAndNTypes)
    return lhs.nsites == rhs.nsites && lhs.Δ == rhs.Δ && typeof(lhs) == typeof(rhs)
end

function Base.:(==)(lhs::SumOfXPhaseTypes, rhs::SumOfXPhaseTypes)
    return lhs.nsites == rhs.nsites && lhs.Ω == rhs.Ω && lhs.ϕ == rhs.ϕ && typeof(lhs) == typeof(rhs)
end

function Base.:(==)(lhs::RydbergHamiltonian, rhs::RydbergHamiltonian)
    return lhs.rydberg_term == rhs.rydberg_term && lhs.rabi_term == rhs.rabi_term && lhs.detuning_term == rhs.detuning_term
end

function Base.:(==)(lhs::RydbergHamiltonian3, rhs::RydbergHamiltonian3)
    return lhs.rydberg_term == rhs.rydberg_term && 
        lhs.rabi_term_hf == rhs.rabi_term_hf && lhs.detuning_term_hf == rhs.detuning_term_hf &&
        lhs.rabi_term_r == rhs.rabi_term_r && lhs.detuning_term_r == rhs.detuning_term_r
end

Base.isreal(::RydInteract) = true
Base.isreal(::SumOfZAndNTypes) = true
Base.isreal(::SumOfXTypes) = true
Base.isreal(::SumOfXPhaseTypes) = false
Base.isreal(h::Add) = all(isreal, subblocks(h))
Base.isreal(h::Scale) = isreal(factor(h)) && isreal(content(h))
Base.isreal(h::RydbergHamiltonian) = !(h.rabi_term isa SumOfXPhase)
Base.isreal(h::RydbergHamiltonian3) = !(h.rabi_term_hf isa SumOfXPhase_01) && !(h.rabi_term_r isa SumOfXPhase_1r)


storage_size(x) = sizeof(x)
function storage_size(h::Hamiltonian)
    return sum(storage_size, h.ts)
end
function storage_size(H::SparseMatrixCSC)
    return sizeof(H.colptr) + sizeof(H.rowval) + sizeof(H.nzval)
end

is_time_dependent(::Nothing) = false
function is_time_dependent(t::Union{RabiTypes, HyperfineRabiTypes, RydbergRabiTypes})
    (t isa Union{SumOfX, SumOfX_01, SumOfX_1r}) && return !is_const_param(t.Ω)
    (t isa Union{SumOfXPhase, SumOfXPhase_01, SumOfXPhase_1r}) && (return !is_const_param(t.Ω) || !is_const_param(t.ϕ))
    return false
end
function is_time_dependent(t::Union{DetuningTypes, HyperfineDetuningTypes, RydbergDetuningTypes})
    (t isa Union{SumOfN, SumOfN_1, SumOfN_r}) && return !is_const_param(t.Δ)
    return false
end
is_time_dependent(h::RydbergHamiltonian) = is_time_dependent(h.rabi_term) || is_time_dependent(h.detuning_term)
is_time_dependent(h::RydbergHamiltonian3) = 
    is_time_dependent(h.rabi_term_hf) || is_time_dependent(h.detuning_term_hf) ||
    is_time_dependent(h.rabi_term_r) || is_time_dependent(h.detuning_term_r)
