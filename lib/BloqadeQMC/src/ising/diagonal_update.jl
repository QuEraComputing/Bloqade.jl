
Base.@propagate_inbounds alignment_check(H::AbstractTFIM, op::NTuple{K, Int}, s1::Bool, s2::Bool) where K =
xor(isferromagnetic(H, getbondsites(H, op)), s1, s2)

Base.@propagate_inbounds alignment_check(H::AbstractLTFIM, op::NTuple{K, Int}, s1::Bool, s2::Bool) where K =
(getoperatortype(H, op) == getbondtype(H, s1, s2))


# insert_diagonal_operator! returns true if operator insertion succeeded
# returns true if operator insertion succeeded
function insert_diagonal_operator!(rng::AbstractRNG, qmc_state::BinaryQMCState{K, V}, H::AbstractIsing{O}, spin_prop::V, n::Int) where {K, V, T, O <: AbstractOperatorSampler{K, T}}
    op = rand(rng, H.op_sampler)
    site1, site2 = getbondsites(H, op)
    @inbounds if issiteoperator(H, op) || alignment_check(H, op, spin_prop[site1], spin_prop[site2])
        qmc_state.operator_list[n] = op
        return op, zero(T)
    else
        return nothing, zero(T)
    end
end

function insert_diagonal_operator!(rng::AbstractRNG, qmc_state::BinaryQMCState{K, V}, H::AbstractLTFIM{<:AbstractImprovedOperatorSampler{K, T}}, spin_prop::V, n::Int) where {K, T, V}
    op = rand(rng, H.op_sampler)

    @inbounds if issiteoperator(H, op)
        qmc_state.operator_list[n] = op
        return op, zero(T)
    else
        t = getoperatortype(H, op)
        lw1 = getlogweight(H, op)
        site1, site2 = getbondsites(H, op)
        real_t = getbondtype(H, spin_prop[site1], spin_prop[site2])
        if t == real_t
            qmc_state.operator_list[n] = op
            return op, lw1
        end

        op2 = convertoperatortype(H, op, real_t)
        lw2 = getlogweight(H, op2)

        if rand(rng) < exp(lw2 - lw1)
            qmc_state.operator_list[n] = op2
            return op2, lw2
        else
            return nothing, zero(lw2)
        end
    end
end

insert_diagonal_operator!(qmc_state, H, spin_prop, n) = insert_diagonal_operator!(Random.GLOBAL_RNG, qmc_state, H, spin_prop, n)
