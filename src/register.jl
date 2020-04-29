export SubspaceReg

struct SubspaceReg{B,ST,SST} <: AbstractRegister{B}
    state::ST
    subspace::SST
    function SubspaceReg{B,ST,SST}(state, subspace) where {B, ST, SST}
        if length(state) != subspace
            DimensionMismatch("size of state $(size(state)) does not match size of subspace $(size(subspace))")
        end
        new{B,ST,SST}(state, subspace)
    end
end

function SubspaceReg(state::VT, subspace::SST) where {VT<:AbstractVector, SST}
    state = reshape(state,:,1)
    SubspaceReg{1,typeof(state),SST}(state, subspace)
end

function SubspaceReg(state::VT, subspace::SST) where {VT<:AbstractMatrix, SST}
    SubspaceReg{size(state,2),VT,SST}(state, subspace)
end

Yao.state(reg::SubspaceReg) = reg.state
Yao.statevec(reg::SubspaceReg) = Yao.matvec(reg.state)
Yao.relaxedvec(reg::SubspaceReg{1}) = vec(reg.state)
Yao.relaxedvec(reg::SubspaceReg) = reg.state
