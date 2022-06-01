"""
    SubspaceArrayReg{D, State, Space} <: AbstractRegister{D}
    SubspaceArrayReg(state, subspace)

Type for registers in a subspace. The subspace must be a
[`Subspace`](@ref).
"""
struct SubspaceArrayReg{D, State<:AbstractVector,Space} <: YaoAPI.AbstractRegister{D}
    natoms::Int
    state::State
    subspace::Space

    function SubspaceArrayReg{D}(state::State, subspace::Subspace) where {State<:AbstractVector, D}
        if length(state) != length(subspace)
            throw(
                DimensionMismatch(
                    "size of state $(size(state, 1)) does not match size of subspace $(length(subspace))",
                ),
            )
        end
        return new{D,State,typeof(subspace)}(subspace.nqubits, state, subspace)
    end
end
SubspaceArrayReg(state::AbstractVector, subspace::Subspace) = SubspaceArrayReg{2}(state, subspace)

Base.copy(reg::SubspaceArrayReg{D}) where D = SubspaceArrayReg{D}(copy(reg.state), copy(reg.subspace))

YaoAPI.nqudits(reg::SubspaceArrayReg) = reg.natoms
YaoAPI.nactive(reg::SubspaceArrayReg) = reg.natoms
YaoArrayRegister.state(reg::SubspaceArrayReg) = reg.state
YaoArrayRegister.statevec(reg::SubspaceArrayReg) = reg.state
YaoArrayRegister.relaxedvec(reg::SubspaceArrayReg) = reg.state
YaoArrayRegister.datatype(reg::SubspaceArrayReg) = eltype(reg.state)

function Adapt.adapt_structure(to, x::SubspaceArrayReg{D}) where D
    return SubspaceArrayReg{D}(Adapt.adapt(to, x.state), x.subspace)
end

"""
    zero_state([T=ComplexF64], n::Int, subspace; nlevel=2)

Create a `SubspaceArrayReg` in zero state in given subspace.

# Arguments

- `T`: optional, element type, default is `ComplexF64`.
- `n`: required, number of atoms (qubits).
- `subspace`: required, the subspace of rydberg state.
"""
YaoArrayRegister.zero_state(subspace::Subspace; nlevel=2) = zero_state(ComplexF64, subspace; nlevel)

function YaoArrayRegister.zero_state(::Type{T}, s::Subspace; nlevel=2) where {T}
    state = zeros(T, length(s))
    state[1] = 1
    return SubspaceArrayReg{nlevel}(state, s)
end

"""
    rand_state([T=ComplexF64], subspace; nlevel=2)

Create a random state in the given subspace.
"""
YaoArrayRegister.rand_state(s::Subspace) = rand_state(ComplexF64, s)

function YaoArrayRegister.rand_state(::Type{T}, s::Subspace; nlevel=2) where {T<:Complex}
    state = normalize!(rand(T, length(s)))
    return SubspaceArrayReg{nlevel}(state, s)
end

"""
    product_state([T=ComplexF64], config, subspace)

Create a product state of given config from `subspace`.
"""
function YaoArrayRegister.product_state(config::DitStr, s::Subspace)
    return YaoArrayRegister.product_state(ComplexF64, config, s)
end

function YaoArrayRegister.product_state(::Type{T}, c::DitStr{D}, s::Subspace) where {T,D}
    c in s.subspace_v || error("$c is not in given subspace")
    state = zeros(T, length(s))
    state[s[c]] = 1
    return SubspaceArrayReg{D}(state, s)
end

# TODO: make upstream implementation more generic
function LinearAlgebra.normalize!(r::SubspaceArrayReg)
    normalize!(r.state)
    return r
end

YaoArrayRegister.isnormalized(r::SubspaceArrayReg) = norm(r.state) ≈ 1
function Base.isapprox(x::SubspaceArrayReg, y::SubspaceArrayReg; kwargs...)
    return nlevel(x) == nlevel(y) && x.natoms == x.natoms && isapprox(x.state, y.state; kwargs...) && (x.subspace == y.subspace)
end

"""
    set_zero_state!(register)

Set the given register to |00...00⟩.
"""
function set_zero_state! end

function set_zero_state!(r::SubspaceArrayReg)
    fill!(r.state, 0)
    r.state[1] = 1
    return r
end

function set_zero_state!(r::ArrayReg)
    fill!(r.state, 0)
    r.state[1, :] .= 1
    return r
end

function Base.:*(bra::YaoArrayRegister.AdjointRegister{D,<:SubspaceArrayReg}, ket::SubspaceArrayReg{D}) where D
    return dot(statevec(parent(bra)), statevec(ket))
end

"""
    space(register)

Return the space of given register.
"""
function space end

space(r::SubspaceArrayReg) = r.subspace
space(r::ArrayReg) = fullspace
space(r::AdjointRegister) = space(parent(r))

# arithmatics operations
# neg
Base.:-(reg::SubspaceArrayReg{D}) where D = SubspaceArrayReg{D}(-reg.state, reg.subspace)

# +, -
for op in [:+, :-]
    @eval function Base.$op(lhs::SubspaceArrayReg{D}, rhs::SubspaceArrayReg{D}) where D
        @assert lhs.natoms == rhs.natoms
        @assert length(lhs.subspace) == length(rhs.subspace)
        return SubspaceArrayReg{D}(($op)(lhs.state, rhs.state), lhs.subspace)
    end
end

function YaoArrayRegister.regadd!(lhs::SubspaceArrayReg{D}, rhs::SubspaceArrayReg{D}) where D
    @assert lhs.natoms == rhs.natoms
    @assert length(lhs.subspace) == length(rhs.subspace)
    lhs.state .+= rhs.state
    return lhs
end

function YaoArrayRegister.regsub!(lhs::SubspaceArrayReg{D}, rhs::SubspaceArrayReg{D}) where D
    @assert lhs.natoms == rhs.natoms
    @assert length(lhs.subspace) == length(rhs.subspace)
    lhs.state .-= rhs.state
    return lhs
end

function YaoArrayRegister.regscale!(lhs::SubspaceArrayReg, x)
    lhs.state .*= x
    return lhs
end

# *, /
for op in [:*, :/]
    @eval function Base.$op(lhs::SubspaceArrayReg{D}, rhs::Number) where D
        return SubspaceArrayReg{D}(($op)(lhs.state, rhs), lhs.subspace)
    end
end

function Base.:*(lhs::Number, rhs::SubspaceArrayReg{D}) where D
    return SubspaceArrayReg{D}(lhs * rhs.state, rhs.subspace)
end

function Base.:(==)(lhs::SubspaceArrayReg, rhs::SubspaceArrayReg)
    return nlevel(lhs) == nlevel(rhs) && lhs.natoms == rhs.natoms && lhs.subspace == rhs.subspace && lhs.state == rhs.state
end

function YaoArrayRegister.most_probable(reg::SubspaceArrayReg{D}, n::Int) where D
    imax = sortperm(abs2.(reg.state); rev = true)[1:n]
    return YaoArrayRegister.DitStr{D, nqubits(reg)}.(reg.subspace.subspace_v[imax])
end
