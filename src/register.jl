export RydbergReg

struct RydbergReg{N,B,ST,SST} <: AbstractRegister{B}
    state::ST
    subspace::SST
    function RydbergReg{N,B,ST,SST}(state::ST, subspace::SST) where {N, B, ST, SST}
        if length(state) != subspace
            DimensionMismatch("size of state $(size(state)) does not match size of subspace $(size(subspace))")
        end
        new{N, B,ST,SST}(state, subspace)
    end
end

function RydbergReg{N}(state::AbstractVector, subspace::SST) where {N, SST}
    state = reshape(state,:,1)
    if !(eltype(subspace) <: BitStr)
        subspace = BitStr64{N}.(subspace)
    end
    return RydbergReg{N, 1, typeof(state), typeof(subspace)}(state, subspace)
end

function RydbergReg{N}(state::VT, subspace::SST) where {N, VT<:AbstractMatrix, SST}
    if !(eltype(subspace) <: BitStr)
        subspace = BitStr64{N}.(subspace)
    end
    return RydbergReg{N, size(state,2),VT,typeof(subspace)}(state, subspace)
end

Yao.nqubits(reg::RydbergReg{N}) where N = N
Yao.nactive(reg::RydbergReg{N}) where N = N
Yao.state(reg::RydbergReg) = reg.state
Yao.statevec(reg::RydbergReg) = Yao.matvec(reg.state)
Yao.relaxedvec(reg::RydbergReg{N, 1}) where N = vec(reg.state)
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
