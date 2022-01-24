"""
    AbstractTerm

Abstract term for hamiltonian terms.
"""
abstract type AbstractTerm end

# TODO: use an actual symbolic zero, maybe consider SymbolicUtils Term
# struct Zero <: Number end

default_unit(unit, x) = x

default_unit(unit, x::Quantity) = uconvert(unit, x).val
default_unit(unit::typeof(NoUnits), x::Quantity) = uconvert(unit, x)
function default_unit(unit, xs::Tuple)
    return map(xs) do x
        default_unit(unit, x)
    end
end

function default_unit(unit, range::AbstractRange)
    a = default_unit(unit, first(range))
    b = default_unit(unit, step(range))
    c = default_unit(unit, last(range))
    return a:b:c
end

function default_unit(unit, x::AbstractArray{S}) where {T, S <: Quantity{T}}
    y = similar(x, T)
    @inbounds for i in eachindex(x)
        y[i] = default_unit(unit, x[i]).val
    end
    return y
end

# function default_unit(unit, f) # parameters are function
#     return function waveform(t)
#         default_unit(unit, f(t))
#     end
# end

const ConstParamType = Union{Number, AbstractVector{<:Number}, NTuple{N, <:Number} where N}
const ConstParamListType = Union{AbstractVector{<:Number}, NTuple{N, <:Number} where N}

assert_has_time_method(::ConstParamType, name) = nothing
assert_has_time_method(::Nothing, name) = nothing # skip symbolic zero (currently nothing)
function assert_has_time_method(fs::Union{AbstractVector, Tuple}, name)
    for f in fs
        assert_has_time_method(f, name)
    end
end

function assert_has_time_method(f, name)
    hasmethod(f, Tuple{Real}) || throw(ArgumentError("invalid input for $name: method $f(::Real) is not defined"))
end

function assert_nsites(nsites::Int, p, name)
    p isa AbstractVector || p isa Tuple || return
    nsites == length(p) ||
        throw(ArgumentError(
            "nsites does not match size of $name " *
            "expect $nsites, got $(length(p))"
    ))
    return
end

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
- `C`: the interaction strength, default unit is `MHz⋅μm^6`. default value is `2π * 858386 * MHz*µm^6`.
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

RydInteract(atoms::AbstractVector) = RydInteract(atoms, 2π * 858386 * MHz*µm^6)

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

# Custom multi-line printing
_print(io::IO, x::AbstractFloat) = printstyled(io, round(x, sigdigits=3); color=:green)
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
    indent = get(io, :indent, 0)
    println(io, " "^indent, nameof(typeof(t)))
    print_term(IOContext(io, :indent=>indent + 1), t)
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
    indent = get(io, :indent, 2)
    for (i, each) in enumerate(t.terms)
        println(io, " "^indent, " Term ", i)
        print_term(IOContext(io, :indent=>indent+2), each)
        
        if i != lastindex(t.terms)
            println(io)
            println(io)
        end
    end
end

function print_term(io::IO, t::Negative)
    print_term(IOContext(io, :negative=>true), t.term)
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
    negative = get(io, :negative, false)
    print(io, " "^indent)
    negative && print(io, "-")
    print(io, "∑(n=1:$nsites) ")
end

function _print_xterm(io::IO, nsites::Int, Ω, ϕ)
    @switch (Ω, ϕ) begin
        @case (::Tuple, ::Tuple) || (::AbstractVector, ::AbstractVector)
            _print_eachterm(io, nsites) do k
                _print_single_xterm(io, Ω[k], ϕ[k])
            end
        @case (_, ::Tuple) || (_, ::AbstractVector)
            _print_eachterm(io, nsites) do k
                _print_single_xterm(io, Ω, ϕ[k])
            end
        @case (::Tuple, _) || (::AbstractVector, _)
            _print_eachterm(io, nsites) do k
                _print_single_xterm(io, Ω[k], ϕ)
            end
        @case _
            _print_sum(io, nsites)
            _print_single_xterm(io, Ω, ϕ)     
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

function SparseArrays.SparseMatrixCSC{Tv, Ti}(term::AbstractTerm, s::AbstractSpace=fullspace) where {Tv, Ti}
    N = space_size(term, s)
    colptr, rowval = sparse_skeleton_csc(Ti, term, s)
    H = SparseMatrixCSC{Tv, Ti}(N, N, colptr, rowval, Vector{Tv}(undef, length(rowval)))
    update_term!(H, term, s)
    return H
end

function SparseArrays.SparseMatrixCSC{Tv, Ti}(term::Hamiltonian, s::AbstractSpace=fullspace) where {Tv, Ti}
    return sum(SparseMatrixCSC{Tv, Ti}(t, s) for t in term.terms)
end

function SparseArrays.SparseMatrixCSC{Tv, Ti}(t::Negative, s::AbstractSpace=fullspace) where {Tv, Ti}
    H = SparseMatrixCSC{Tv, Ti}(t.term, s)
    lmul!(-one(Tv), H.nzval)
    return H
end

# NOTE: we use BlasInt as Ti, since MKLSparse uses BlasInt as Ti
# on some devices BlasInt != Int, thus this is necessary to trigger
# dispatch to MKLSparse
SparseArrays.SparseMatrixCSC{Tv}(term::AbstractTerm, s::AbstractSpace=fullspace) where Tv = SparseMatrixCSC{Tv, BlasInt}(term, s)
SparseArrays.SparseMatrixCSC(term::AbstractTerm, s::AbstractSpace=fullspace) = SparseMatrixCSC{Complex{eltype(term)}}(term, s)

sparse_skeleton_csc(t::AbstractTerm, s::AbstractSpace=fullspace) = sparse_skeleton_csc(Int, t, s)

function sparse_skeleton_csc(::Type{Ti}, t::XTerm, s::FullSpace) where Ti
    natoms = nsites(t)
    sz = hilbert_space(natoms)
    colptr = Ti[natoms*i+1 for i in 0:length(sz)]
    rowval = Vector{Ti}(undef, natoms * length(sz))
    @inbounds @batch for col in sz
        colrowval = map(1:natoms) do k
            rhs = flip(col, 1 << (k - 1))
            return rhs + 1
        end
        rowval[natoms*col+1:natoms*col+natoms] = sort!(colrowval)
    end
    return colptr, rowval
end

function sparse_skeleton_csc(::Type{Ti}, t::XTerm, s::Subspace) where Ti
    colptr = Vector{Ti}(undef, length(s)+1)
    rowval_list = Vector{Vector{Ti}}(undef, length(s))
    colptr[1] = 1

    @inbounds @batch for i in 1:length(s)
        lhs = s.subspace_v[i]
        colrowval = Ti[]
        for k in 1:nsites(t)
            rhs = flip(lhs, 1 << (k - 1))
            if haskey(s, rhs)
                push!(colrowval, s[rhs])
            end
        end
        sort!(colrowval)
        rowval_list[i] = colrowval
    end

    rowval = Vector{Ti}(undef, sum(length, rowval_list))
    ptr = 1
    for col in 1:length(rowval_list) # O(length(s))
        rowval_col = rowval_list[col]
        rowval_nnz = length(rowval_col)
        rowval[ptr:ptr+rowval_nnz-1] = rowval_col
        colptr[col+1] = colptr[col] + rowval_nnz
        ptr += rowval_nnz
    end
    return colptr, rowval
end

# Diagonal
function sparse_skeleton_csc(::Type{Ti}, t::Union{RydInteract, ZTerm, NTerm}, s::AbstractSpace) where Ti
    m = space_size(t, s)
    colptr = Vector{Ti}(1:(m+1))
    rowval = Vector{Ti}(1:m)
    return colptr, rowval
end

"""
    update_term!(H, term[, space=fullspace])

Update matrix `H` based on the given Hamiltonian term. This can be faster when the sparse structure of
`H` is known (e.g `H` is a `SparseMatrixCSC`).
"""
function update_term! end

Base.@propagate_inbounds function foreach_nnz(f, H::SparseMatrixCSC)
    @batch for lhs in 1:size(H, 1)
        @inbounds start = H.colptr[lhs]
        @inbounds stop = H.colptr[lhs+1]-1

        for k in start:stop
            @inbounds rhs = H.rowval[k]
            f(k, lhs, rhs)
        end
    end
end

function update_term!(H::AbstractSparseMatrix, t::AbstractTerm, ::FullSpace=fullspace)
    nzval = nonzeros(H)
    foreach_nnz(H) do k, col, row
        @inbounds nzval[k] = term_value(t, col-1, row-1, col, row)
    end
    return H
end

function update_term!(H::AbstractSparseMatrix, t::AbstractTerm, s::Subspace)
    nzval = nonzeros(H)
    @inbounds foreach_nnz(H) do k, col, row
        lhs = s.subspace_v[col]
        rhs = s.subspace_v[row]
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

Base.@propagate_inbounds function term_value(t::OrNegative{<:XTerm}, lhs, rhs, col, row)
    col == row && return zero(eltype(t))
    mask = lhs ⊻ rhs
    l = unsafe_log2i(mask) + 1
    l_site = rhs & mask
    return getterm(t, l, l_site)
end

Base.@propagate_inbounds function term_value(t::OrNegative{<:Union{ZTerm, NTerm}}, lhs, rhs, col, row)
    col != row && return zero(eltype(t))
    sigma_z = zero(eltype(t))
    for i in 1:nsites(t)
        sigma_z += getterm(t, i, readbit(lhs, i))
    end
    return sigma_z
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


function rydberg_h(atoms, C, Ω, ϕ, Δ)
    return RydInteract(atoms, C) + XTerm(length(atoms), Ω, ϕ) - NTerm(length(atoms), Δ)
end

function rydberg_h(atoms, Ω, ϕ, Δ)
    return RydInteract(atoms) + XTerm(length(atoms), Ω, ϕ) - NTerm(length(atoms), Δ)
end

"""
    rydberg_h(atoms; [C=2π * 858386 * MHz*µm^6], Ω[, ϕ, Δ])

Create a rydberg hamiltonian

```math
∑ \\frac{C}{|r_i - r_j|^6} n_i n_j + \\frac{Ω}{2} σ_x - Δ σ_n
```

shorthand for

```julia
RydInteract(C, atoms) + XTerm(length(atoms), Ω, ϕ) + ZTerm(length(atoms), Δ)
```

# Arguments

- `atoms`: a collection of atom positions

# Keyword Arguments

- `C`: optional, default unit is `MHz*µm^6`, interation parameter,
    see also [`RydInteract`](@ref).
- `Ω`: required, default unit is `MHz`, Rabi frequencies, see [`XTerm`](@ref).
- `Δ`: optional, default unit is `MHz`, detuning parameter, see [`NTerm`](@ref).
- `ϕ`: optional, does not have unit, the phase, see [`XTerm`](@ref).

!!! tips

    The parameters of Hamiltonian have their own default units to match hardware,
    one can use [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)
    to specify their units explicitly. If the units are specified explicitly,
    they will be converted to default units automatically.

# Example

```julia-repl
julia> using EaRyd

julia> atoms = generate_sites(SquareLattice(), 3, 3);

julia> rydberg_h(atoms; Δ=1.2, Ω=1.1)
Hamiltonian
  Term 1
   ∑(n=1:9) 1.1/2 σ^x

  Term 2
   -∑(n=1:9) 1.2 n

  Term 3
   ∑(n=1:9) 686.0/|r_i - r_j|^6 n_i n_j
```

```julia-repl
julia> rydberg_h(atoms; Δ=1.2, Ω=1.1, ϕ=2.1)
Hamiltonian
  Term 1
   ∑(n=1:9) 1.1/2 (e^{2.1i}|0)⟨1| + e^{-2.1i}|1⟩⟨0|)

  Term 2
   -∑(n=1:9) 1.2 n

  Term 3
   ∑(n=1:9) 686.0/|r_i - r_j|^6 n_i n_j
```
"""
function rydberg_h(atoms; Ω, Δ=nothing, C=nothing, ϕ=nothing) 
    if eltype(atoms) <: Tuple
        atoms = map(atoms) do each
            RydAtom(each)
        end
        atoms = SVector{length(atoms)}(atoms)
    end

    if ϕ === nothing
        term = XTerm(length(atoms), Ω)
    else
        term = XTerm(length(atoms), Ω, ϕ)
    end

    if Δ !== nothing
        term -= NTerm(length(atoms), Δ)
    end

    if C === nothing
        term += RydInteract(atoms)
    else
        term += RydInteract(atoms, C)
    end
    return term
end

function is_time_dependent(h::XTerm)
    return !(h.Ωs isa ConstParamType || isnothing(h.Ωs)) ||
        !(h.ϕs isa ConstParamType || isnothing(h.ϕs))
end

function is_time_dependent(h::Union{ZTerm, NTerm})
    return !(h.Δs isa ConstParamType)
end

# NOTE: we currently assume atom positions are constant
# might change in the future

is_time_dependent(h::RydInteract) = false
is_time_dependent(t::Negative) = is_time_dependent(t.term)

"""
    is_time_dependent(h::Hamiltonian)

Check if a hamiltonian is time-dependent.
"""
function is_time_dependent(h::Hamiltonian)
    any(is_time_dependent, h.terms)
end

# Time dependent Term
attime(x::Nothing, t::Real) = nothing
attime(x::Number, t::Real) = x
attime(x, t::Real) = x(t)
attime(x::AbstractArray, t::Real) = attime.(x, t)
attime(x::Tuple, t::Real) = attime.(x, t)

function (tm::XTerm)(t::Real)
    return XTerm(tm.nsites, attime(tm.Ωs, t), attime(tm.ϕs, t))
end

function (tm::ZTerm)(t::Real)
    return ZTerm(tm.nsites, attime(tm.Δs, t))
end

function (tm::NTerm)(t::Real)
    return NTerm(tm.nsites, attime(tm.Δs, t))
end

(tm::Negative)(t::Real) = Negative(tm.term(t))

function (tm::Hamiltonian)(t::Real)
    return Hamiltonian(map(x->x(t), tm.terms))
end

# fallback to constants
(tm::AbstractTerm)(t::Real) = tm

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
