distance(a, b) = sqrt(mapreduce(x -> x^2, +, a .- b))

function YaoAPI.mat(::Type{T}, h::AbstractBlock, space::FullSpace) where {T}
    return YaoAPI.mat(T, h)
end

YaoAPI.mat(::Type{T}, h::XPhase{P}) where {T, P} = PermMatrix([2, 1], T[exp(h.ϕ * im), exp(-h.ϕ * im)])
YaoAPI.mat(::Type{T}, h::XPhase_01{P}) where {T, P} = PermMatrix([2, 1, 3], T[exp(h.ϕ * im), exp(-h.ϕ * im), 0])
YaoAPI.mat(::Type{T}, h::XPhase_1r{P}) where {T, P} = PermMatrix([1, 3, 2], T[0, exp(h.ϕ * im), exp(-h.ϕ * im)])

YaoAPI.mat(::Type{T}, h::PuPhase{P}) where {T, P} = sparse([1], [2], T[exp(h.ϕ * im)], 2, 2)
YaoAPI.mat(::Type{T}, h::PuPhase_1r{P}) where {T, P} = sparse([2], [3], T[exp(h.ϕ * im)], 3, 3)
YaoAPI.mat(::Type{T}, h::PuPhase_01{P}) where {T, P} = sparse([1], [2], T[exp(h.ϕ * im)], 3, 3)

YaoAPI.mat(::Type{T}, h::PdPhase{P}) where {T, P} = sparse([2], [1], T[exp(-h.ϕ * im)], 2, 2)
YaoAPI.mat(::Type{T}, h::PdPhase_01{P}) where {T, P} = sparse([2], [1], T[exp(-h.ϕ * im)], 3, 3)
YaoAPI.mat(::Type{T}, h::PdPhase_1r{P}) where {T, P} = sparse([3], [2], T[exp(-h.ϕ * im)], 3, 3)

YaoAPI.mat(::Type{T}, ::OpX_01) where T = PermMatrix([2, 1, 3], T[1, 1, 0])
YaoAPI.mat(::Type{T}, ::OpX_1r) where T = PermMatrix([1, 3, 2], T[0, 1, 1])

YaoAPI.mat(::Type{T}, ::OpZ_01) where T = Diagonal(T[1, -1, 0])
YaoAPI.mat(::Type{T}, ::OpZ_1r) where T = Diagonal(T[0, 1, -1])

YaoAPI.mat(::Type{T}, ::OpN_1) where T = sparse([2], [2], T[1], 3, 3)
YaoAPI.mat(::Type{T}, ::OpN_r) where T = sparse([3], [3], T[1], 3, 3)

YaoAPI.mat(::Type{T}, ::OpPu_01) where T = sparse([1], [2], T[1], 3, 3)
YaoAPI.mat(::Type{T}, ::OpPu_1r) where T = sparse([2], [3], T[1], 3, 3)

YaoAPI.mat(::Type{T}, ::OpPd_01) where T = sparse([2], [1], T[1], 3, 3)
YaoAPI.mat(::Type{T}, ::OpPd_1r) where T = sparse([3], [2], T[1], 3, 3)

function YaoAPI.mat(::Type{T}, h::AbstractTerm{D}, space::FullSpace = fullspace) where {T, D}
    blocks = YaoBlocks.Optimise.simplify(h; rules = [YaoBlocks.Optimise.to_basictypes])
    return YaoAPI.mat(T, blocks)
end

function YaoAPI.mat(::Type{T}, h::AbstractTerm{D}, space::Subspace) where {T, D}
    blocks = YaoBlocks.Optimise.simplify(h; rules = [YaoBlocks.Optimise.to_basictypes])
    return YaoAPI.mat(T, blocks, space)
end

# fallback
YaoBlocks.mat(::Type{T}, x::AbstractBlock, s::Subspace) where {T} = _toregular_matrix(mat(T, x))[s, s]
_toregular_matrix(x::AbstractMatrix) = x
_toregular_matrix(x::Union{PermMatrix,IMatrix}) = SparseMatrixCSC(x)
YaoBlocks.mat(x::AbstractBlock, s::Subspace) = mat(promote_type(ComplexF64, YaoBlocks.parameters_eltype(x)), x, s)

YaoBlocks.mat(::Type{T}, d::TrivialGate{N}, s::Subspace) where {T,N} = IMatrix{T}(length(s))
YaoBlocks.mat(::Type{T}, pb::PutBlock{N,C}, s::Subspace) where {T,N,C} =
    _cunmat(s.subspace_v, s.map, (), (), mat(T, pb.content), pb.locs)
YaoBlocks.mat(::Type{T}, c::Subroutine{D,<:AbstractBlock}, s::Subspace) where {D,T} =
    mat(T, PutBlock(c.n, c.content, c.locs), s)
YaoBlocks.mat(::Type{T}, x::Scale, s::Subspace) where {T} = YaoBlocks.factor(x) * mat(T, content(x), s)
YaoBlocks.mat(::Type{T}, blk::Daggered, s::Subspace) where {T} = adjoint(mat(T, content(blk), s))
YaoBlocks.mat(::Type{T}, c::ControlBlock{N,BT,C}, s::Subspace) where {T,N,BT,C} =
    _cunmat(s.subspace_v, s.map, c.ctrl_locs, c.ctrl_config, mat(T, c.content), c.locs)
YaoBlocks.mat(::Type{T}, x::Add{N}, s::Subspace) where {N,T} = mapreduce(x -> mat(T, x, s), +, x.list)
function YaoBlocks.mat(::Type{T}, c::ChainBlock, s::Subspace) where {T}
    # @warn "ChainBlock is applying subspace on its subblocks, which may be inconsistent with multiply and then take subspace."
    return prod(x -> mat(T, x, s), Iterators.reverse(c.blocks))
end
function YaoBlocks.mat(::Type{T}, te::TimeEvolution{N}, s::Subspace) where {T,N}
    @warn "TimeEvolution block is applying subspace on its subblock, which may be inconsistent with time evolve and then take subspace."
    return exp(-im * T(te.dt) * Matrix(mat(T, te.H, s)))
end
# NOTE: CachedBlock is not yet implemented.

function _cunmat(
    subspace_v::AbstractVector{TI},
    map,
    cbits::NTuple{C,Int},
    cvals::NTuple{C,Int},
    U0::AbstractMatrix{T},
    locs::NTuple{M,Int},
) where {TI,T,C,M}
    U = staticize(all(diff(collect(locs)) .> 0) ? U0 : reorder(U0, collect(locs) |> sortperm))
    # ctest = YaoBlocks.controller(TI, cbits, cvals)
    # NOTE: remove the following lines after the update of BitBasis!
    do_mask = bmask(TI, cbits...)
    target =
        length(cvals) == 0 ? zero(TI) :
        mapreduce(xy -> (xy[2] == 1 ? one(TI) << TI(xy[1] - 1) : zero(TI)), |, zip(cbits, cvals))
    ctest = b -> ismatch(b, do_mask, target)
    return _main_loop(subspace_v, map, U, ctest, locs)
end

function _main_loop(subspace_v, map, U::AbstractMatrix{T}, ctest, locs) where {T}
    N = length(subspace_v)
    MM = size(U, 1)
    colptr = Vector{Int}(undef, N + 1)
    colptr[1] = 1
    locs_zeroclear_mask = ~bmask(locs...)
    NNZ_MAX = max(1, maximum([count(!iszero, view(U, :, j)) for j in 1:size(U, 2)])) * N
    rowval = Vector{Int}(undef, NNZ_MAX)
    nzval = Vector{T}(undef, NNZ_MAX)
    nel = 0

    @inbounds for j in 1:N
        J = subspace_v[j]
        jindex = readbit(J, locs...)
        if ctest(J)
            I0 = J & locs_zeroclear_mask
            nj = 0
            for iindex in 0:MM-1
                mij = U[iindex+1, jindex+1]
                iszero(mij) && continue
                I = setlocbits(I0, locs, iindex)
                i = get(map, I, -1)
                if i > 0
                    nel += 1
                    rowval[nel] = i
                    nzval[nel] = mij
                    nj += 1
                end
            end
            colptr[j+1] = colptr[j] + nj
        else
            nel += 1
            colptr[j+1] = colptr[j] + 1
            rowval[nel] = j
            nzval[nel] = one(T)
        end
    end
    return SparseMatrixCSC(N, N, colptr, rowval[1:nel], nzval[1:nel])
end

@inline @generated function setlocbits(I0, locs::NTuple{N}, iindex) where {N}
    quote
        @nexprs $N i -> I0 = I0 | (readbit(iindex, i) << (locs[i] - 1))
        return I0
    end
end
