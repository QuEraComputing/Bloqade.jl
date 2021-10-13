using Yao.YaoArrayRegister, Yao.YaoBase
using Base.Cartesian: @nexprs

function YaoBase.instruct!(
    r::RydbergReg{ComplexLayout,VT},
    op,
    locs::Tuple,
    control_locs::Tuple,
    control_bits::Tuple,
) where {T,VT<:AbstractArray{T}}
    mask = bmask(control_locs)
    onemask = bmask(Int[control_locs[i] for i=1:length(control_locs) if control_bits[i]==1])
    op = YaoArrayRegister.sort_unitary(_mat(T, op), locs)
    res = zero(r.state)
    mp = r.subspace.map   # this is slower than regular dict!
    _instruct!(res, r.state, r.subspace.subspace_v, mp, YaoArrayRegister.staticize(op), locs, mask, onemask)
    r.state .= res
    return r
end

function YaoBase.instruct!(
    r::RydbergReg{ComplexLayout,VT},
    op,
    locs::Tuple,
    ::Tuple{},
    ::Tuple{},
) where {T,VT<:AbstractArray{T}}
    op = YaoArrayRegister.sort_unitary(_mat(T, op), locs)
    res = zero(r.state)
    mp = r.subspace.map   # this is slower than regular dict!
    _instruct!(res, r.state, r.subspace.subspace_v, mp, YaoArrayRegister.staticize(op), locs)
    r.state .= res
    return r
end

for OP in YaoArrayRegister.SPECIALIZATION_LIST
    @eval _mat(::Type{T}, ::Val{$(QuoteNode(OP))}) where T = mat(T, Yao.ConstGate.$OP)
end
_mat(::Type{T}, m::AbstractMatrix{T}) where T = m

function YaoBase.instruct!(r::RydbergReg, op, locs::Tuple)
    instruct!(r, op, locs, (), ())
end

function YaoBase.instruct!(
    r::RydbergReg{L,VT},
    op,
    locs::Tuple,
    control_locs::Tuple,
    control_bits::Tuple,
    theta::Number,
) where {L,T,VT<:AbstractArray{T}}
    m = YaoArrayRegister.rot_mat(T, op, theta)
    instruct!(r, m, locs, control_locs, control_bits)
end

function YaoBase.instruct!(r::RydbergReg, op, locs::Tuple, theta::Number)
    instruct!(r, op, locs, (), (), theta)
end

function _instruct!(res, state, subspace_v, map, m, locs)
    locs_zeroclear_mask = ~bmask(locs...)
    @inbounds for j=1:length(subspace_v)
        J = subspace_v[j]
        jindex = readbit(J, locs...)
        # zero clear bits at locs
        I0 = J & locs_zeroclear_mask
        for iindex in 0:size(m,1)-1
            mij = m[iindex+1,jindex+1]
            iszero(mij) && continue
            I = setlocbits(I0, locs, iindex)
            i = get(map, I, -1)
            if i > 0
                res[i] += mij * state[j]
            end
        end
    end
    return res
end

function _instruct!(res, state, subspace_v, map, m, locs, mask, onemask)
    locs_zeroclear_mask = ~bmask(locs...)
    @inbounds for j=1:length(subspace_v)
        J = subspace_v[j]
        if testval(J, mask, onemask)
            jindex = readbit(J, locs...)
            # zero clear bits at locs
            I0 = J & locs_zeroclear_mask
            for iindex in 0:size(m,1)-1
                mij = m[iindex+1,jindex+1]
                iszero(mij) && continue
                I = setlocbits(I0, locs, iindex)
                i = get(map, I, -1)
                if i > 0
                    res[i] += mij * state[j]
                end
            end
        else
            res[j] = state[j]
        end
    end
    return res
end

@inline @generated function setlocbits(I0, locs::NTuple{N}, iindex) where N
    quote
        @nexprs $N i->I0 = I0 | (readbit(iindex, i) << (locs[i]-1))
        return I0
    end
end
