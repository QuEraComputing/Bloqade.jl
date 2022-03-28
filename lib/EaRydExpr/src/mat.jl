distance(a, b) = sqrt(mapreduce(x->x^2, +, a .- b))

function YaoAPI.mat(::Type{T}, h::AbstractBlock, space::FullSpace) where T
    return YaoAPI.mat(T, h)
end

function YaoAPI.mat(::Type{T}, h::XPhase) where T
    return PermMatrix([2, 1], [exp(h.ϕ*im), exp(-h.ϕ*im)])
end

function YaoAPI.mat(::Type{T}, h::AbstractTerm, space::FullSpace=fullspace) where T
    return YaoAPI.mat(T, emit_lowered(h))
end

function YaoAPI.mat(::Type{T}, h::AbstractTerm, space::Subspace) where T
    return YaoAPI.mat(T, emit_lowered(h), space)
end

# fallback
YaoBlocks.mat(::Type{T}, x::AbstractBlock, s::Subspace) where T =
    _toregular_matrix(mat(T, x))[s, s]
_toregular_matrix(x::AbstractMatrix) = x
_toregular_matrix(x::Union{PermMatrix, IMatrix}) = SparseMatrixCSC(x)
YaoBlocks.mat(x::AbstractBlock, s::Subspace) =
    mat(promote_type(ComplexF64, YaoBlocks.parameters_eltype(x)), x, s)

YaoBlocks.mat(::Type{T}, d::TrivialGate{N}, s::Subspace) where {T,N} = IMatrix{length(s),T}()
YaoBlocks.mat(::Type{T}, pb::PutBlock{N,C}, s::Subspace) where {T,N,C} =
    _cunmat(s.subspace_v, s.map, (), (), mat(T, pb.content), pb.locs)
YaoBlocks.mat(::Type{T}, c::Subroutine{D, <:AbstractBlock}, s::Subspace) where {D,T} =
    mat(T, PutBlock(c.n, c.content, c.locs), s)
YaoBlocks.mat(::Type{T}, x::Scale, s::Subspace) where {T} = YaoBlocks.factor(x) * mat(T, content(x), s)
YaoBlocks.mat(::Type{T}, blk::Daggered, s::Subspace) where {T} = adjoint(mat(T, content(blk), s))
YaoBlocks.mat(::Type{T}, c::ControlBlock{N,BT,C}, s::Subspace) where {T,N,BT,C} =
    _cunmat(s.subspace_v, s.map, c.ctrl_locs, c.ctrl_config, mat(T, c.content), c.locs)
YaoBlocks.mat(::Type{T}, x::Add{N}, s::Subspace) where {N,T} =
    mapreduce(x -> mat(T, x, s), +, x.list)
function YaoBlocks.mat(::Type{T}, c::ChainBlock, s::Subspace) where {T}
    @warn "ChainBlock is applying subspace on its subblocks, which may be inconsistent with multiply and then take subspace."
    prod(x -> mat(T, x, s), Iterators.reverse(c.blocks))
end
function YaoBlocks.mat(::Type{T}, te::TimeEvolution{N}, s::Subspace) where {T,N}
    @warn "TimeEvolution block is applying subspace on its subblock, which may be inconsistent with time evolve and then take subspace."
    exp(-im * T(te.dt) * Matrix(mat(T, te.H, s)))
end
# NOTE: CachedBlock is not yet implemented.

function _cunmat(
    subspace_v,
    map,
    cbits::NTuple{C,Int},
    cvals::NTuple{C,Int},
    U0::AbstractMatrix{T},
    locs::NTuple{M,Int},
) where {T,C,M}
    U = staticize(all(diff(collect(locs)) .> 0) ? U0 : reorder(U0, collect(locs) |> sortperm))
    ctest = YaoBlocks.controller(cbits, cvals)
    return _main_loop(subspace_v, map, U, ctest, locs)
end

function _main_loop(subspace_v, map, U::AbstractMatrix{T}, ctest, locs) where T
    N = length(subspace_v)
    MM = size(U, 1)
    colptr = Vector{Int}(undef, N + 1); colptr[1] = 1
    locs_zeroclear_mask = ~bmask(locs...)
    NNZ_MAX = max(1, maximum([count(!iszero, view(U, :, j)) for j=1:size(U,2)])) * N
    rowval = Vector{Int}(undef, NNZ_MAX)
    nzval = Vector{T}(undef, NNZ_MAX)
    nel = 0

    @inbounds for j=1:N
        J = subspace_v[j]
        jindex = readbit(J, locs...)
        if ctest(J)
            I0 = J & locs_zeroclear_mask
            nj = 0
            for iindex in 0:MM-1
                mij = U[iindex+1,jindex+1]
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

@inline @generated function setlocbits(I0, locs::NTuple{N}, iindex) where N
    quote
        @nexprs $N i->I0 = I0 | (readbit(iindex, i) << (locs[i]-1))
        return I0
    end
end


# function YaoAPI.mat(::Type{T}, h::RydInteract) where T
#     n = nqubits(h)
#     values = Vector{T}(undef, 2^n)
#     for idx in 1:2^n
#         v = zero(T)
#         for i in 1:n, j in 1:i-1
#             if (readbit(idx, i) == 1) && (readbit(idx, j) == 1)
#                 r_i, r_j = h.atoms[i], h.atoms[j]
#                 alpha = h.C / distance(r_i, r_j)^6
#                 v += alpha
#             end
#         end
#         values[idx] = v
#     end
#     return Diagonal(values)
# end

# function YaoAPI.mat(::Type{T}, h::SumOfX) where T
#     is_time_function(h.Ω) && error("cannot get matrix of a time-dependent operator")
#     ex = if h.Ω isa Vector
#         sum(h.Ω[i] * put(h.nsites, i=>X) for i in 1:h.nsites)
#     else
#         h.Ω * sum(put(h.nsites, i=>X) for i in 1:h.nsites)
#     end
#     return mat(T, ex)
# end

# function YaoAPI.mat(::Type{T}, h::Union{SumOfN, SumOfZ}) where T
#     is_time_function(h.Δ) && error("cannot get matrix of a time-dependent operator")
#     op = h isa SumOfN ? ConstGate.P1 : Z
#     ex = if h.Δ isa Vector
#         sum(h.Δ[i] * put(h.nsites, i=>op) for i in 1:h.nsites)
#     else
#         h.Δ * sum(put(h.nsites, i=>op) for i in 1:h.nsites)
#     end
#     return mat(T, ex)
# end

# function YaoAPI.mat(::Type{T}, h::SumOfXPhase) where T
#     (is_time_function(h.Ω) || is_time_function(h.ϕ)) &&
#         error("cannot get matrix of a time-dependent operator")

#     ex = @switch (h.Ω, h.ϕ) begin
#         @case (::Vector, ::Vector)
#             sum(h.Ω[i] * put(h.nsites, i=>XPhase(h.ϕ[i])) for i in 1:h.nsites)
#         @case (Ω::Vector, ϕ)
#             sum(Ω_i * put(h.nsites, i=>XPhase(ϕ)) for (i, Ω_i) in enumerate(Ω))
#         @case (Ω, ϕ::Vector)
#             sum(Ω * put(h.nsites, i=>XPhase(ϕ_i)) for (i, ϕ_i) in enumerate(ϕ))
#         @case (Ω, ϕ)
#             sum(Ω, put(h.nsites, i=>XPhase(ϕ)) for i in 1:h.nsites)
#     end

#     return mat(T, ex)
# end
