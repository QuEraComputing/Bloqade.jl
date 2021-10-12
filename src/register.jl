abstract type MemoryLayout end
struct RealLayout <: MemoryLayout end
struct ComplexLayout <: MemoryLayout end

struct RydbergReg{Layout <: MemoryLayout, State, Space <: Subspace} <: Yao.AbstractRegister{1}
    natoms::Int
    layout::Layout
    state::State
    subspace::Space

    function RydbergReg{Layout, State, Space}(natoms, layout, state, subspace) where {Layout, State, Space}
        new{Layout, State, Space}(natoms, layout, state, subspace)
    end

    function RydbergReg(natoms::Int, layout::MemoryLayout, state::State, subspace::Subspace) where {State <: AbstractArray}
        if size(state, 1) != length(subspace)
            throw(DimensionMismatch("size of state $(size(state, 1)) does not match size of subspace $(length(subspace))"))
        end
        new{typeof(layout), State, typeof(subspace)}(natoms, layout, state, subspace)
    end
end

function Adapt.adapt_structure(to, x::RydbergReg)
    return RydbergReg(x.natoms, x.layout, adapt(to, x.state), x.subspace)
end

MemoryLayout(::Yao.ArrayReg) = ComplexLayout()
MemoryLayout(::RydbergReg{Layout}) where Layout = Layout()

Base.copy(reg::RydbergReg{L, State, Space}) where {L, State, Space} =
    RydbergReg{L, State, Space}(reg.natoms, reg.layout, copy(reg.state), copy(reg.subspace))

RydbergReg{Layout}(reg::RydbergReg{Layout}) where {Layout <: MemoryLayout} = copy(reg)

function RydbergReg{RealLayout}(reg::RydbergReg{ComplexLayout})
    T = real(eltype(reg.state))
    dst = similar(reg.state, T, (length(reg.state), 2))
    copyto!(dst, reinterpret(reshape, T, reg.state)')
    return RydbergReg(reg.natoms, RealLayout(), dst, reg.subspace)
end

function RydbergReg{ComplexLayout}(reg::RydbergReg{RealLayout})
    @views state = reg.state[:, 1] + reg.state[:, 2] * im
    return RydbergReg(reg.natoms, ComplexLayout(), state, reg.subspace)
end

"""
    RydbergReg(natoms::Int, state::AbstractVector, subspace::Subspace)

Create a `RydbergReg` from state vector and its corresponding subspace of `natoms`.
"""
function RydbergReg(natoms::Int, state::AbstractVector, subspace::Subspace) where {N}
    return RydbergReg(natoms, ComplexLayout(), reshape(state, size(state, 1)), subspace)
end

"""
    RydbergReg(natoms::Int, state::AbstractMatrix, subspace::Subspace)

Create a `RydbergReg` from real value storage.
"""
function RydbergReg(natoms::Int, state::AbstractMatrix{<:Real}, subspace::Subspace)
    return RydbergReg(natoms, RealLayout(), state, subspace)
end

function RydbergReg(natoms::Int, state::AbstractMatrix{<:Complex}, subspace::Subspace)
    error("support of batched register is dropped, please open an issue for your use case")
end

Yao.nqubits(reg::RydbergReg) = reg.natoms
Yao.nactive(reg::RydbergReg) = reg.natoms
Yao.state(reg::RydbergReg) = reg.state
Yao.statevec(reg::RydbergReg) = reg.state
Yao.relaxedvec(reg::RydbergReg) = reg.state
Yao.datatype(reg::RydbergReg{ComplexLayout}) = eltype(reg.state)
Yao.datatype(reg::RydbergReg{RealLayout}) = Complex{eltype(reg.state)}

"""
    zero_state([T=ComplexF64], n::Int, subspace[, layout=ComplexLayout()])

Create a `RydbergReg` in zero state in given subspace.

# Arguments

- `T`: optional, element type, default is `ComplexF64`.
- `n`: required, number of atoms (qubits).
- `subspace`: required, the subspace of rydberg state.
- `layout`: optional, memory layout, default is `ComplexLayout`.

# Memory Layout

When the hamiltonian is a real hermitian, it can be more efficient
to use the `RealLayout` over `ComplexLayout` which stores the complex
-value state vector as as `length(state)×2` matrix, the first column
is the real component and the second column is the imaginary component.
"""
Yao.zero_state(n::Int, subspace::Subspace, layout=ComplexLayout()) = zero_state(ComplexF64, n, subspace, layout)
Yao.zero_state(::Type{T}, natoms::Int, s::Subspace) where T = zero_state(T, natoms, s, ComplexLayout())

function Yao.zero_state(::Type{T}, natoms::Int, s::Subspace, layout::ComplexLayout) where {T <: Complex}
    state = zeros(T, length(s))
    state[1] = 1
    return RydbergReg(natoms, layout, state, s)
end

function Yao.zero_state(::Type{T}, natoms::Int, s::Subspace, layout::RealLayout) where {T <: Complex}
    state = zeros(real(T), length(s), 2)
    state[1, 1] = 1
    return RydbergReg(natoms, layout, state, s)
end

function Yao.rand_state(natoms::Int, s::Subspace, layout::MemoryLayout = ComplexLayout())
    Yao.rand_state(ComplexF64, natoms, s, layout)
end

function Yao.rand_state(::Type{T}, natoms::Int, s::Subspace, layout::ComplexLayout) where {T <: Complex}
    state = normalize!(rand(T, length(s)))
    return RydbergReg(natoms, layout, state, s)
end

function Yao.rand_state(::Type{T}, natoms::Int, s::Subspace, layout::RealLayout) where {T <: Complex}
    state = normalize!(rand(real(T), length(s), 2))
    return RydbergReg(natoms, layout, state, s)
end

function Yao.product_state(natoms::Int, config::BitStr, s::Subspace, layout::MemoryLayout=ComplexLayout())
    Yao.product_state(ComplexF64, natoms, config, s, layout)
end

function Yao.product_state(::Type{T}, natoms::Int, c::BitStr, s::Subspace, layout::ComplexLayout) where T
    c in s.subspace_v || error("$c is not in given subspace")
    state = zeros(T, length(s))
    state[s[c]] = 1
    return RydbergReg(natoms, layout, state, s)
end

function Yao.product_state(::Type{T}, natoms::Int, c::BitStr, s::Subspace, layout::RealLayout) where T
    c in s.subspace_v || error("$c is not in given subspace")
    state = zeros(T, length(s), 2)
    state[s[c], 1] = 1
    return RydbergReg(natoms, layout, state, s)
end

# TODO: make upstream implementation more generic
function LinearAlgebra.normalize!(r::RydbergReg)
    normalize!(r.state)
    return r
end

Yao.isnormalized(r::RydbergReg) = norm(r.state) ≈ 1
function Base.isapprox(x::RydbergReg{L}, y::RydbergReg{L}; kwargs...) where L
    return isapprox(x.state, y.state; kwargs...) && (x.subspace == y.subspace)
end
function Base.isapprox(x::RydbergReg{RealLayout}, y::RydbergReg{ComplexLayout}; kw...)
    return isapprox(x.state[:, 1] + im * x.state[:, 2], y.state; kw...)
end
function Base.isapprox(x::RydbergReg{ComplexLayout}, y::RydbergReg{RealLayout}; kw...)
    return isapprox(y, x; kw...)
end

"""
    set_zero_state!(register)

Set the given register to |00...00⟩.
"""
function set_zero_state! end

function set_zero_state!(r::RydbergReg)
    fill!(r.state, 0)
    r.state[1] = 1
    return r
end

function set_zero_state!(r::Yao.ArrayReg)
    fill!(r.state, 0)
    r.state[1, :] .= 1
    return r
end

function Base.:*(bra::Yao.YaoArrayRegister.AdjointRegister{1, <:RydbergReg}, ket::RydbergReg)
    return dot(statevec(parent(bra)), statevec(ket))
end
