# Todo: Get rid of "improved" prefix

# abstract type AbstractOperatorSampler{K, T, P <: AbstractProbabilityVector{T}} end            We probably don't need Abstract type
firstindex(::OperatorSampler) = 1               # changed argument from taking abstract type to concrete type
lastindex(os::OperatorSampler) = length(os)
@inline normalization(os::OperatorSampler) = normalization(os.pvec)

# abstract type AbstractOperatorSampler{K, T, P} <: AbstractOperatorSampler{K, T, P} end        Again, probably no need for Abstract type

struct OperatorSampler{K, T, P} # <: OperatorSampler{K, T, P}                   no abstract subtyping
    operators::Vector{NTuple{K, Int}}
    pvec::ProbabilityAlias              # replaced P with ProbabilityAlias 
    op_log_weights::Vector{T}
end

# What do I do with OperatorSampler in argument here?
function OperatorSampler(H::Type{<:Hamiltonian{2, <:OperatorSampler}}, operators::Vector{NTuple{K, Int}}, p::Vector{T}) where {T <: AbstractFloat, K}
    @assert length(operators) == length(p) "Given vectors must have the same length!"

    op_log_weights = log.(p)

    max_mel_ops = Vector{NTuple{K, Int}}()
    p_modified = Vector{T}()

    # fill with all the site operators first
    for (i, op) in enumerate(operators)
        if issiteoperator(H, op) && isdiagonal(H, op)
            push!(max_mel_ops, op)
            push!(p_modified, @inbounds p[i])
        end
    end
    idx = findall(isbondoperator(H), operators)
    ops = operators[idx]
    p_mod = p[idx]

    perm = sortperm(ops, by=getbondsites(H))
    ops = ops[perm]
    p_mod = p_mod[perm]

    op_groups = Dict{NTuple{2, Int}, Vector{NTuple{K, Int}}}()
    p_groups = Dict{NTuple{2, Int}, Vector{T}}()

    while !isempty(ops)
        op = pop!(ops)
        site1, site2 = getbondsites(H, op)
        p_ = pop!(p_mod)

        if haskey(op_groups, (site1, site2))
            push!(op_groups[(site1, site2)], op)
            push!(p_groups[(site1, site2)], p_)
        else
            op_groups[(site1, site2)] = [op]
            p_groups[(site1, site2)] = [p_]
        end
    end

    for (k, p_gr) in p_groups
        i = argmax(p_gr)
        push!(max_mel_ops, op_groups[k][i])
        push!(p_modified, p_gr[i])
    end

    pvec = probability_vector(p_modified)
    return OperatorSampler{K, T, typeof(pvec)}(max_mel_ops, pvec, op_log_weights)                   # Here, we return only the max matrix elements and associated compressed probability vector. Full information still contained in op_log_weights
end

@inline rand(rng::AbstractRNG, os::OperatorSampler) = @inbounds os.operators[rand(rng, os.pvec)]

Base.@propagate_inbounds getlogweight(os::OperatorSampler, w::Int) = os.op_log_weights[w]

@inline length(os::OperatorSampler) = length(os.operators)

##############################################################################
