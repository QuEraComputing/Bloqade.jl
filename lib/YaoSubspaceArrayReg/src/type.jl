"""
    SubspaceArrayReg <: AbstractRegister{2}
    SubspaceArrayReg(state, subspace)

Type for registers in a subspace. The subspace must be a
[`Subspace`](@ref).
"""
struct SubspaceArrayReg{State <: AbstractVector, Space} <: YaoAPI.AbstractRegister{2}
    natoms::Int
    state::State
    subspace::Space

    function SubspaceArrayReg(state::State, subspace::Subspace) where {State <: AbstractVector}
        if length(state) != length(subspace)
            throw(DimensionMismatch("size of state $(size(state, 1)) does not match size of subspace $(length(subspace))"))
        end
        new{State, typeof(subspace)}(subspace.nqubits, state, subspace)
    end
end

Base.copy(reg::SubspaceArrayReg) = SubspaceArrayReg(copy(reg.state), copy(reg.subspace))

YaoAPI.nqudits(reg::SubspaceArrayReg) = reg.natoms
YaoAPI.nactive(reg::SubspaceArrayReg) = reg.natoms
YaoArrayRegister.state(reg::SubspaceArrayReg) = reg.state
YaoArrayRegister.statevec(reg::SubspaceArrayReg) = reg.state
YaoArrayRegister.relaxedvec(reg::SubspaceArrayReg) = reg.state
YaoArrayRegister.datatype(reg::SubspaceArrayReg) = eltype(reg.state)

"""
    zero_state([T=ComplexF64], n::Int, subspace)

Create a `SubspaceArrayReg` in zero state in given subspace.

# Arguments

- `T`: optional, element type, default is `ComplexF64`.
- `n`: required, number of atoms (qubits).
- `subspace`: required, the subspace of rydberg state.
"""
YaoArrayRegister.zero_state(subspace::Subspace) = zero_state(ComplexF64, subspace)

function YaoArrayRegister.zero_state(::Type{T}, s::Subspace) where T
    state = zeros(T, length(s))
    state[1] = 1
    return SubspaceArrayReg(state, s)
end

"""
    rand_state(subspace)

Create a random state in the given subspace.
"""
YaoArrayRegister.rand_state(s::Subspace) = rand_state(ComplexF64, s)

function YaoArrayRegister.rand_state(::Type{T}, s::Subspace) where {T <: Complex}
    state = normalize!(rand(T, length(s)))
    return SubspaceArrayReg(state, s)
end

"""
    product_state(config, subspace)

Create a product state of given config from `subspace`.
"""
function YaoArrayRegister.product_state(config::BitStr, s::Subspace)
    YaoArrayRegister.product_state(ComplexF64, config, s)
end

function YaoArrayRegister.product_state(::Type{T}, c::BitStr, s::Subspace) where T
    c in s.subspace_v || error("$c is not in given subspace")
    state = zeros(T, length(s))
    state[s[c]] = 1
    return SubspaceArrayReg(state, s)
end

# TODO: make upstream implementation more generic
function LinearAlgebra.normalize!(r::SubspaceArrayReg)
    normalize!(r.state)
    return r
end

YaoArrayRegister.isnormalized(r::SubspaceArrayReg) = norm(r.state) ≈ 1
function Base.isapprox(x::SubspaceArrayReg, y::SubspaceArrayReg; kwargs...)
    return x.natoms == x.natoms && isapprox(x.state, y.state; kwargs...) && (x.subspace == y.subspace)
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

function Base.:*(bra::YaoArrayRegister.AdjointRegister{2, <:SubspaceArrayReg}, ket::SubspaceArrayReg)
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
Base.:-(reg::SubspaceArrayReg) = SubspaceArrayReg(-reg.state, reg.subspace)

# +, -
for op in [:+, :-]
    @eval function Base.$op(lhs::SubspaceArrayReg, rhs::SubspaceArrayReg)
        @assert lhs.natoms == rhs.natoms
        @assert length(lhs.subspace) == length(rhs.subspace)
        return SubspaceArrayReg(($op)(lhs.state, rhs.state), lhs.subspace)
    end
end

function YaoArrayRegister.regadd!(lhs::SubspaceArrayReg, rhs::SubspaceArrayReg)
    @assert lhs.natoms == rhs.natoms
    @assert length(lhs.subspace) == length(rhs.subspace)
    lhs.state .+= rhs.state
    return lhs
end

function YaoArrayRegister.regsub!(lhs::SubspaceArrayReg, rhs::SubspaceArrayReg)
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
    @eval function Base.$op(lhs::SubspaceArrayReg, rhs::Number)
        return SubspaceArrayReg(($op)(lhs.state, rhs), lhs.subspace)
    end
end

function Base.:*(lhs::Number, rhs::SubspaceArrayReg)
    return SubspaceArrayReg(lhs * rhs.state, rhs.subspace)
end

function Base.:(==)(lhs::SubspaceArrayReg, rhs::SubspaceArrayReg)
    lhs.natoms == rhs.natoms &&
    lhs.subspace == rhs.subspace &&
    lhs.state == rhs.state
end