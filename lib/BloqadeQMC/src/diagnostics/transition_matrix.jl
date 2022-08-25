using OnlineStats


abstract type AbstractTransitionMatrix{T}; end

# returns the transition matrix with row sums normalized to 1
function LinearAlgebra.normalize(M::AbstractTransitionMatrix{T}) where T
    rowsums = sum(M.m, dims=2)

    normalized = similar(M.m, float(T))
    for i in eachindex(rowsums)
        if !iszero(rowsums[i])
            @. normalized[i, :] = M.m[i, :] / rowsums[i]
        else
            @. normalized[i, :] = float(M.m[i, :])
        end
    end
    return normalized
end


struct NoTransitionMatrix <: AbstractTransitionMatrix{Int}; end
@inline OnlineStats.fit!(M::NoTransitionMatrix, _, _) = M
LinearAlgebra.normalize(M::NoTransitionMatrix) = nothing

# For now we'll assume that our operators are at most 2-local


# This transition matrix is constructed by assuming that the tuple (k, abs(t), w)
#  uniquely defines the matrix element; later on this may no longer be true
#  (i.e. for very large systems) since we may need to start mapping multiple
#  different operator tuples acting on different sites to the same weight index.
struct PositionDependentTransitionMatrix{T <: Integer} <: AbstractTransitionMatrix{T}
    m::Matrix{T}
    nspins::Int

    function PositionDependentTransitionMatrix{T}(nspins::Int, w_range::UnitRange{Int}) where T
        m = zeros(T, nspins + last(w_range), nspins + last(w_range))
        new{T}(m, nspins)
    end
end
PositionDependentTransitionMatrix(nspins, w_range) =
    PositionDependentTransitionMatrix{Int}(nspins, w_range)
@inline function OnlineStats.fit!(M::PositionDependentTransitionMatrix, init::NTuple{K, Int}, final::NTuple{K, Int}) where K
    ns = M.nspins
    if init[1] == 1
        i = 2*init[3] - (init[2] > 0)
        f = 2*final[3] - (final[2] > 0)
    else init[2] == 2
        i = init[3] + ns
        f = final[3] + ns
    end
    M.m[f, i] += 1
    return M
end

# This transition matrix only tracks transitions between different matrix element
#  types, ignoring the index w, it is therefore position independent.
struct TransitionMatrix{T <: Integer} <: AbstractTransitionMatrix{T}
    m::Matrix{T}

    function TransitionMatrix{T}() where T
        m = zeros(T, 6, 6)
        new{T}(m)
    end
end
TransitionMatrix() = TransitionMatrix{Int}()

@inline function OnlineStats.fit!(M::TransitionMatrix, init::NTuple{K, Int}, final::NTuple{K, Int}) where K
    if init[1] == 1
        M.m[abs(final[2]), abs(init[2])] += 1
    else
        M.m[final[2] + 2, init[2] + 2] += 1
    end
    return M
end


struct CombinedTransitionMatrix{T <: Integer} <: AbstractTransitionMatrix{T}
    pdmatrix::PositionDependentTransitionMatrix{T}
    tmatrix::TransitionMatrix{T}

    function CombinedTransitionMatrix{T}(nspins::Int, w_range::UnitRange{Int}) where T
        new{T}(PositionDependentTransitionMatrix{T}(nspins, w_range), TransitionMatrix{T}())
    end
end
CombinedTransitionMatrix(nspins, w_range) = CombinedTransitionMatrix{Int}(nspins, w_range)
@inline function OnlineStats.fit!(M::CombinedTransitionMatrix, init::NTuple{K, Int}, final::NTuple{K, Int}) where K
    fit!(M.pdmatrix, init, final)
    fit!(M.tmatrix, init, final)
    return M
end
