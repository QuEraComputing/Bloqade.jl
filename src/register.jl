struct RydbergReg{N,B,ST,SST} <: Yao.AbstractRegister{B}
    state::ST
    subspace::SST
    function RydbergReg{N,B,ST,SST}(state::ST, subspace::SST) where {N, B, ST, SST}
        if size(state, 1) != length(subspace)
            throw(DimensionMismatch("size of state $(size(state, 1)) does not match size of subspace $(length(subspace))"))
        end
        new{N, B,ST,SST}(state, subspace)
    end
end

function RydbergReg{N}(state::AbstractVector, subspace::Subspace) where {N}
    state = reshape(state,:,1)
    return RydbergReg{N, 1, typeof(state), typeof(subspace)}(state, subspace)
end

function RydbergReg{N}(state::VT, subspace::Subspace) where {N, VT<:AbstractMatrix}
    return RydbergReg{N, size(state,2),VT,typeof(subspace)}(state, subspace)
end

Yao.nqubits(::RydbergReg{N}) where N = N
Yao.nactive(::RydbergReg{N}) where N = N
Yao.state(reg::RydbergReg) = reg.state
Yao.statevec(reg::RydbergReg) = Yao.YaoArrayRegister.matvec(reg.state)
Yao.relaxedvec(reg::RydbergReg{N, 1}) where N = vec(reg.state)
Yao.relaxedvec(reg::RydbergReg) = reg.state

Base.copy(reg::RydbergReg{N, B, ST, SST}) where {N, B, ST, SST} =
    RydbergReg{N, B, ST, SST}(copy(reg.state), copy(reg.subspace))

"""
    zero_state([T=ComplexF64], n::Int, subspace; nbatch=1)

Create a `RydbergReg` in zero state in given subspace.
"""
Yao.zero_state(n::Int, subspace::Subspace; nbatch=1) = zero_state(ComplexF64, n, subspace; nbatch=nbatch)

function Yao.zero_state(::Type{T}, n::Int, s::Subspace; nbatch=1) where {T <: Complex}
    st = zeros(T, length(s), nbatch)
    st[1, :] .= 1
    return RydbergReg{n}(st, s)
end

Yao.rand_state(n::Int, s::Subspace; nbatch::Int=1) = Yao.rand_state(ComplexF64, n, s; nbatch=nbatch)
function Yao.rand_state(::Type{T}, n::Int, s::Subspace; nbatch::Int=1) where T
    st = Yao.batch_normalize!(rand(T, length(s), nbatch))
    return RydbergReg{n}(st, s)
end

Yao.product_state(n::Int, c::BitStr, s::Subspace; nbatch::Int=1) = Yao.product_state(ComplexF64, n, c, s; nbatch=nbatch)
function Yao.product_state(::Type{T}, n::Int, c::BitStr, s::Subspace; nbatch::Int=1) where T
    c in s.subspace_v || error("$c is not in given subspace")
    st = zeros(T, length(s), nbatch)
    st[s[c], :] .= 1
    return RydbergReg{n}(st, s)
end

# TODO: make upstream implementation more generic
Yao.isnormalized(r::RydbergReg) = all(sum(abs2, r.state, dims = 1) .≈ 1)

Base.isapprox(x::RydbergReg, y::RydbergReg; kwargs...) = isapprox(x.state, y.state; kwargs...) && (x.subspace == y.subspace)

"""
    set_zero_state!(register)

Set the given register to |00...00⟩.
"""
function set_zero_state! end

set_zero_state!(r::RydbergReg) = (_set_zero_state!(r.state); r)
set_zero_state!(r::Yao.ArrayReg) = (_set_zero_state!(r.state); r)

function _set_zero_state!(st::AbstractMatrix{T}) where T
    fill!(st, zero(T))
    st[1, :] .= 1
    return st
end
