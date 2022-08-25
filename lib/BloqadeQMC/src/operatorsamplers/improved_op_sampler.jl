abstract type AbstractImprovedOperatorSampler{K, T, P} <: AbstractOperatorSampler{K, T, P} end

struct ImprovedOperatorSampler{K, T, P} <: AbstractImprovedOperatorSampler{K, T, P}
    operators::Vector{NTuple{K, Int}}
    pvec::P
    op_log_weights::Vector{T}
end

# only supports the LTFIM/Rydberg cases for now
function ImprovedOperatorSampler(H::Type{<:Hamiltonian{2, <:AbstractOperatorSampler}}, operators::Vector{NTuple{K, Int}}, p::Vector{T}) where {T <: AbstractFloat, K}
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
    return ImprovedOperatorSampler{K, T, typeof(pvec)}(max_mel_ops, pvec, op_log_weights)
end

@inline rand(rng::AbstractRNG, os::ImprovedOperatorSampler) = @inbounds os.operators[rand(rng, os.pvec)]

Base.@propagate_inbounds getlogweight(os::ImprovedOperatorSampler, w::Int) = os.op_log_weights[w]

@inline length(os::ImprovedOperatorSampler) = length(os.operators)

##############################################################################
