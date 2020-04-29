export RydbergReg

struct RydbergReg{N,B,ST,SST} <: AbstractRegister{B}
    state::ST
    subspace::SST
    function RydbergReg{N,B,ST,SST}(state, subspace) where {N, B, ST, SST}
        if length(state) != subspace
            DimensionMismatch("size of state $(size(state)) does not match size of subspace $(size(subspace))")
        end
        new{N, B,ST,SST}(state, subspace)
    end
end

function RydbergReg{N}(state::VT, subspace::SST) where {N, VT<:AbstractVector, SST}
    state = reshape(state,:,1)
    return RydbergReg{N, 1, VT, SST}(state, subspace)
end

function RydbergReg{N}(state::VT, subspace::SST) where {N, VT<:AbstractMatrix, SST}
    return RydbergReg{N, size(state,2),VT,SST}(state, subspace)
end

Yao.nqubits(reg::RydbergReg{N}) where N = N
Yao.nactive(reg::RydbergReg{N}) where N = N
Yao.state(reg::RydbergReg) = reg.state
Yao.statevec(reg::RydbergReg) = Yao.matvec(reg.state)
Yao.relaxedvec(reg::RydbergReg{1}) = vec(reg.state)
Yao.relaxedvec(reg::RydbergReg) = reg.state

"""
    zero_state([T=ComplexF64], n::Int, subspace; nbatch=1)

Create a Rydberg zero state in given subspace.
"""
zero_state(n::Int, subspace; nbatch=1) = zero_state(ComplexF64, n, subspace; nbatch=nbatch)

function zero_state(::Type{T}, n::Int, subspace; nbatch=1) where T
    st = zeros(T, length(subspace), nbatch)
    st[1, :] .= 1
    return RydbergReg{n}(st, subspace)
end

# TODO: make upstream implementation more generic
Yao.isnormalized(r::RydbergReg) = all(sum(abs2, r.state, dims = 1) .≈ 1)
