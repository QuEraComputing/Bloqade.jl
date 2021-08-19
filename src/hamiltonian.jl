"""
    AbstractTerm

Abstract term for hamiltonian terms.
"""
abstract type AbstractTerm end

default_unit(unit, x) = x

default_unit(unit, x::Quantity) = uconvert(unit, x).val
default_unit(unit::typeof(NoUnits), x::Quantity) = uconvert(unit, x)
default_unit(unit, xs::Tuple{}) = ()
default_unit(unit, xs::Tuple{T}) where T = default_unit(unit, first(xs))
default_unit(unit, xs::Tuple) = (default_unit(unit, Base.heads(xs)), default_unit(unit, Base.tail(xs))...)

function default_unit(unit, x::AbstractArray{S}) where {T, S <: Quantity{T}}
    y = similar(x, T)
    @inbounds for i in eachindex(x)
        y[i] = default_unit(unit, x[i]).val
    end
    return y
end

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

    function XTerm{Omega, Phi}(nsites::Int, Ωs::Omega, ϕs::Phi) where {Omega, Phi}
        new{Omega, Phi}(nsites, Ωs, ϕs)
    end

    function XTerm(nsites::Int, Ωs, ϕs)
        Ωs = default_unit(MHz, Ωs)
        ϕs = default_unit(NoUnits, ϕs)
        new{typeof(Ωs), typeof(ϕs)}(nsites, Ωs, ϕs)
    end
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

    function ZTerm{Delta}(nsites::Int, Δs::Delta) where {Delta}
        new{Delta}(nsites, Δs)
    end

    function ZTerm(nsites::Int, Δs)
        Δs = default_unit(MHz, Δs)
        new{typeof(Δs)}(nsites, Δs)
    end
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

    function NTerm{Delta}(nsites::Int, Δs::Delta) where {Delta}
        new{Delta}(nsites, Δs)
    end

    function NTerm(nsites::Int, Δs)
        Δs = default_unit(MHz, Δs)
        new{typeof(Δs)}(nsites, Δs)
    end
end

struct Hamiltonian{Terms <: Tuple} <: AbstractTerm
    terms::Terms
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
nsites(t::RydInteract) = length(t.atoms)

hilbert_space(n::Int) = 0:((1<<n)-1)
hilbert_space(t::AbstractTerm) = hilbert_space(nsites(t))

Base.eltype(t::XTerm) = eltype(t.Ωs)
Base.eltype(t::ZTerm) = eltype(t.Δs)
Base.eltype(t::NTerm) = eltype(t.Δs)
Base.eltype(t::RydInteract) = typeof(t.C)
Base.eltype(t::Hamiltonian) = promote_type(eltype.(t.terms)...)

Base.getindex(t::Hamiltonian, i::Int) = t.terms[i]

# Custom multi-line printing
_print(io::IO, x::AbstractFloat) = printstyled(io, @sprintf("%.2f", x); color=:green)
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
    println(io, nameof(typeof(t)))
    print_term(IOContext(io, :indent=>1), t)
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
    for (i, each) in enumerate(t.terms)
        println(io, " Term ", i)
        print_term(IOContext(io, :indent=>2), each)
        
        if i != lastindex(t.terms)
            println(io)
            println(io)
        end
    end
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
    print(io, " "^indent, "∑(n=1:$nsites) ")
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

Base.:(+)(x::AbstractTerm, y::AbstractTerm) = Hamiltonian((x, y))
Base.:(+)(x::AbstractTerm, y::Hamiltonian) = Hamiltonian((x, y.terms...))
Base.:(+)(x::Hamiltonian, y::AbstractTerm) = Hamiltonian((x.terms..., y))
Base.:(+)(x::Hamiltonian, y::Hamiltonian) = Hamiltonian((x.terms..., y.terms...))

# absorb - to RHS
Base.:(-)(x::AbstractTerm, y::AbstractTerm) = Hamiltonian((x, -y))
Base.:(-)(x::AbstractTerm, y::Hamiltonian) = Hamiltonian((x, map(-, y.terms)...))
Base.:(-)(x::Hamiltonian, y::AbstractTerm) = Hamiltonian((x.terms..., -y))
Base.:(-)(x::Hamiltonian, y::Hamiltonian) = Hamiltonian((x.terms..., map(-, y.terms)...))

Base.:(-)(x::XTerm{<:ConstParamType}) = XTerm(x.nsites, map(-, x.Ωs), x.ϕs)
Base.:(-)(x::XTerm) = XTerm(x.nsites, t->-x.Ωs(t), x.ϕs)

Base.:(-)(x::ZTerm{<:ConstParamType}) = ZTerm(x.nsites, map(-, x.Δs))
Base.:(-)(x::ZTerm) = ZTerm(x.nsites, t->-x.Δs(t))
Base.:(-)(x::NTerm{<:ConstParamType}) = NTerm(x.nsites, map(-, x.Δs))
Base.:(-)(x::NTerm) = NTerm(x.nsites, t->-x.Δs(t))


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

function getterm(t::XTerm, k, k_site)
    if k_site == 0
        return getscalarmaybe(t.Ωs, k)/2 * exp(im * getscalarmaybe(t.ϕs, k))
    else
        return getscalarmaybe(t.Ωs, k)/2 * exp(-im * getscalarmaybe(t.ϕs, k))
    end
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

"""
    to_matrix!(dst, term[, subspace])

Create given term to a matrix and assign it to `dst`. An optional argument `subspace`
can be taken to construct the matrix in subspace. `dst` should be initialized to zero
entries by the user.
"""
function to_matrix! end

function SparseArrays.SparseMatrixCSC{Tv}(term::AbstractTerm) where {Tv}
    N = 1 << nsites(term)
    H = SparseMatrixCOO{Tv}(undef, N, N)
    to_matrix!(H, term)
    return SparseMatrixCSC(H)
end

function SparseArrays.SparseMatrixCSC{Tv}(term::AbstractTerm, s::Subspace) where {Tv}
    N = length(s)
    H = SparseMatrixCOO{Tv}(undef, N, N)
    to_matrix!(H, term, s)
    return SparseMatrixCSC(H)
end

SparseArrays.SparseMatrixCSC(term::AbstractTerm, xs...) = SparseMatrixCSC{ComplexF64}(term, xs...)

Yao.mat(::Type{T}, t::AbstractTerm) where T = SparseMatrixCSC{T}(t)
Yao.mat(::Type{T}, t::AbstractTerm, s::Subspace) where T = SparseMatrixCSC{T}(t, s)
Yao.mat(t::AbstractTerm, xs...) = Yao.mat(ComplexF64, t, xs...)

# full space
# C/|r_i - r_j|^6 n_i n_j
# specialize on COO
function to_matrix!(dst::SparseMatrixCOO, t::RydInteract)
    n = nsites(t)
    @inbounds for i in 1:n, j in 1:i-1
        r_i, r_j = t.atoms[i], t.atoms[j]
        alpha = t.C / distance(r_i, r_j)^6
        # |11⟩⟨11|
        for lhs in itercontrol(nsites(t), [i, j], [1, 1])
            dst[lhs+1, lhs+1] = alpha # this will be accumulated by COO format converter
        end
    end
    return dst
end

# Ω ⋅ (e^{iϕ}|0)⟨1| + e^{-iϕ} |1⟩⟨0|)
function to_matrix!(dst::AbstractMatrix{T}, t::XTerm) where T
    @inbounds for lhs in hilbert_space(t)
        for k in 1:nsites(t)
            k_site = readbit(lhs, k)
            rhs = flip(lhs, 1 << (k - 1))
            dst[lhs+1, rhs+1] = getterm(t, k, k_site)
        end
    end
    return dst
end

function to_matrix!(dst::AbstractMatrix{T}, t::ZTerm) where T
    @inbounds for lhs in hilbert_space(t)
        sigma_z = zero(T)
        for k in 1:nsites(t)
            sigma_z += getterm(t, k, readbit(lhs, k))
        end
        dst[lhs+1, lhs+1] = sigma_z
    end
    return dst
end

function to_matrix!(dst::AbstractMatrix{T}, t::NTerm) where T
    @inbounds for lhs in hilbert_space(t)
        n_ij = zero(T)
        for k in 1:nsites(t)
            n_ij += getterm(t, k, readbit(lhs, k))
        end
        dst[lhs+1, lhs+1] = n_ij
    end
    return dst
end

function to_matrix!(dst::AbstractMatrix{T}, t::Hamiltonian, xs...) where T
    for term in t.terms
        to_matrix!(dst, term, xs...)
    end
    return dst
end

# subspace
function to_matrix!(dst::AbstractMatrix, t::XTerm, s::Subspace)
    @inbounds for (lhs, i) in s
        for k in 1:nsites(t)
            k_site = readbit(lhs, k)
            rhs = flip(lhs, 1 << (k - 1))
            if haskey(s, rhs)
                dst[i, s[rhs]] = getterm(t, k, k_site)
            end
        end
    end
    return dst
end

function to_matrix!(dst::AbstractMatrix{T}, t::ZTerm, s::Subspace) where T
    @inbounds for (lhs, i) in s
        sigma_z = zero(T)
        for k in 1:nsites(t)
            sigma_z += getterm(t, k, readbit(lhs, k))
        end
        dst[i, i] = sigma_z
    end
    return dst
end

function to_matrix!(dst::AbstractMatrix{T}, t::NTerm, s::Subspace) where T
    @inbounds for (lhs, i) in s
        n_ij = zero(T)
        for k in 1:nsites(t)
            n_ij += getterm(t, k, readbit(lhs, k))
        end
        dst[i, i] = n_ij
    end
    return dst
end


function to_matrix!(dst::SparseMatrixCOO, t::RydInteract, s::Subspace)
    n = nsites(t)
    for (k, lhs) in enumerate(s.subspace_v)
        for i in 1:n, j in 1:i-1
            if (readbit(lhs, i) == 1) && (readbit(lhs, j) == 1)
                r_i, r_j = t.atoms[i], t.atoms[j]
                alpha = t.C / distance(r_i, r_j)^6
                dst[k, k] = alpha
            end
        end
    end
    return dst
end


"""
    update_term!(H, term[, subspace])

Update matrix `H` based on the given Hamiltonian term. This can be faster when the sparse structure of
`H` is known (e.g `H` is a `SparseMatrixCSC`). It fallbacks to `to_matrix!(H, term[, subspace])` if the
sparse structure is unknown.
"""
function update_term! end

# forward to to_matrix! as fallback
update_term!(H::AbstractMatrix, t::AbstractTerm, s::Nothing) = update_term!(H, t)
update_term!(H::SparseMatrixCSC, t::AbstractTerm, s::Nothing) = update_term!(H, t) # disambiguity
update_term!(H::AbstractMatrix, t::AbstractTerm) = to_matrix!(H, t)
update_term!(H::AbstractMatrix, t::AbstractTerm, s::Subspace) = to_matrix!(H, t, s)

# specialize on sparse matrix
update_term!(H::AbstractSparseMatrix, t::AbstractTerm, s::Subspace) = update_term!(H, t, s.subspace_v)

function foreach_nnz(f, H::SparseMatrixCSC)
    for lhs in 1:size(H, 1)
        @inbounds start = H.colptr[lhs]
        @inbounds stop = H.colptr[lhs+1]-1

        for k in start:stop
            @inbounds rhs = H.rowval[k]
            f(k, lhs, rhs)
        end
    end
end

function update_term!(H::AbstractSparseMatrix, t::AbstractTerm)
    nzval = nonzeros(H)
    foreach_nnz(H) do k, col, row
        @inbounds nzval[k] = term_value(t, col-1, row-1, col, row)
    end
    return H
end

function update_term!(H::AbstractSparseMatrix, t::AbstractTerm, subspace_v::AbstractVector)
    nzval = nonzeros(H)
    @inbounds foreach_nnz(H) do k, col, row
        lhs = subspace_v[col]
        rhs = subspace_v[row]

        nzval[k] = term_value(t, lhs, rhs, col, row)
    end
    return H
end

"""
    term_value(term, lhs, rhs, col, row)

Return the value of given term at `H[col, row]` with left basis `lhs` and right basis `rhs`.
For full space, `lhs = col - 1` and `rhs = row - 1`, for subspace, `lhs = subspace_v[col]` and
`rhs = subspace_v[row]`.
"""
function term_value end

@generated function term_value(t::Hamiltonian{Term}, lhs, rhs, col, row) where Term
    ex = Expr(:block)

    push!(ex.args, Expr(:meta, :inline, :propagate_inbounds))
    push!(ex.args, :(val = term_value(t.terms[1], lhs, rhs, col, row)))
    for k in 2:length(Term.parameters)
        push!(ex.args, :(val += term_value(t.terms[$k], lhs, rhs, col, row)))
    end

    push!(ex.args, :val)
    return ex
end

Base.@propagate_inbounds function term_value(t::XTerm, lhs, rhs, col, row)
    col == row && return zero(eltype(t))
    mask = lhs ⊻ rhs
    l = unsafe_log2i(mask) + 1
    l_site = rhs & mask
    return getterm(t, l, l_site)
end

Base.@propagate_inbounds function term_value(t::ZTerm, lhs, rhs, col, row)
    col != row && return zero(eltype(t))
    sigma_z = zero(eltype(t))
    for i in 1:nsites(t)
        sigma_z += getterm(t, i, readbit(lhs, i))
    end
    return sigma_z
end

Base.@propagate_inbounds function term_value(t::NTerm, lhs, rhs, col, row)
    col != row && return zero(eltype(t))
    n_ij = zero(eltype(t))
    for i in 1:nsites(t)
        n_ij += getterm(t, i, readbit(lhs, i))
    end
    return n_ij
end


Base.@propagate_inbounds function term_value(t::RydInteract, lhs, rhs, col, row)
    col != row && return zero(eltype(t))
    # all the nonzeros indices contains
    v = zero(eltype(t))
    n = nsites(t)
    for i in 1:n, j in 1:i-1
        if (readbit(lhs, i) == 1) && (readbit(lhs, j) == 1)
            r_i, r_j = t.atoms[i], t.atoms[j]
            alpha = t.C / distance(r_i, r_j)^6
            v += alpha
        end
    end
    return v
end

Base.@propagate_inbounds getscalarmaybe(x::AbstractVector, k) = x[k]
Base.@propagate_inbounds getscalarmaybe(x::Tuple, k) = x[k]
@inline getscalarmaybe(x::Number, k) = x
@inline getscalarmaybe(x::Nothing, k) = 0
@inline getscalarmaybe(x, k) = x

"""
    simple_rydberg(n::Int, ϕ::Number)

Create a simple rydberg hamiltonian that has only [`XTerm`](@ref).
"""
simple_rydberg(n::Int, ϕ::Number) = XTerm(n, one(ϕ), ϕ)

"""
    rydberg_h(C, atoms, Ω, ϕ, Δ)

Create a rydberg hamiltonian, shorthand for
`RydInteract(C, atoms) + XTerm(length(atoms), Ω, ϕ) + ZTerm(length(atoms), Δ)`

```math
∑ \\frac{C}{|r_i - r_j|^6} n_i n_j + Ω σ_x + Δ σ_z
```
"""
function rydberg_h(atoms, C, Ω, ϕ, Δ)
    return RydInteract(atoms, C) + XTerm(length(atoms), Ω, ϕ) + ZTerm(length(atoms), Δ)
end

# Time dependent Term
attime(x::Nothing, t::Real) = nothing
attime(x::Number, t::Real) = x
attime(x, t::Real) = x(t)
attime(x::AbstractArray, t::Real) = attime.(x, t)

function (tm::XTerm)(t::Real)
    return XTerm(tm.nsites, attime(tm.Ωs, t), attime(tm.ϕs, t))
end

function (tm::ZTerm)(t::Real)
    return ZTerm(tm.nsites, attime(tm.Δs, t))
end

function (tm::NTerm)(t::Real)
    return NTerm(tm.nsites, attime(tm.Δs, t))
end

function (tm::Hamiltonian)(t::Real)
    return Hamiltonian(map(x->x(t), tm.terms))
end

# fallback to constants
(tm::AbstractTerm)(t::Real) = tm
