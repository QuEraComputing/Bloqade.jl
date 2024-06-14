@inline function multibranch_kernel!(qmc_state::BinaryQMCState, H::AbstractIsing, ccount::Int, leg::Int, a::Int)
    Ns = nspins(H)
    LegType, Associates, leg_sites = qmc_state.leg_types, qmc_state.associates, qmc_state.leg_sites
    in_cluster, cstack, current_cluster = qmc_state.in_cluster, qmc_state.cstack, qmc_state.current_cluster

    @inbounds ll, la = LegType[leg], LegType[a]
    @inbounds sl, sa = leg_sites[leg], leg_sites[a]
    # TODO: check if this is inputting the spins in the correct order
    if sl > sa
        preflip_bond_type = getbondtype(H, ll, la)
        postflip_bond_type = getbondtype(H, !ll, !la)
    else
        preflip_bond_type = getbondtype(H, la, ll)
        postflip_bond_type = getbondtype(H, !la, !ll)
    end

    # now check all associates and add them to the cluster
    @inbounds while a != 0 && iszero(in_cluster[a])
        push!(cstack, a)
        in_cluster[a] = ccount
        push!(current_cluster, a)
        a = Associates[a]
    end

    return preflip_bond_type, postflip_bond_type
end

multibranch_acceptance(H::AbstractIsing, lnA::T) where {T <: Real} =
    haslongitudinalfield(H) ? exp(min(lnA, zero(lnA))) : T(0.5)
multibranch_acceptance(H::AbstractRydberg, lnA::T) where {T <: Real} = exp(min(lnA, zero(T)))
# in the TFIM case, acceptance rate is exactly 1 so we set it to 1/2 to ensure ergodicity
multibranch_acceptance(H::AbstractTFIM, lnA::T) where {T <: Real} = T(0.5)


function multibranch_cluster_update!(rng::AbstractRNG, lsize::Int, qmc_state::BinaryQMCState, H::AbstractIsing, d::Diagnostics)
    return generic_cluster_update!(rng, multibranch_kernel!, multibranch_acceptance, lsize, qmc_state, H, d) # call cluster_update with the proper kernel
end

# Here, I cut the simplified implementation for TFIM 

#############################################################################


function multibranch_update!(rng::AbstractRNG, qmc_state::BinaryQMCState, H::AbstractIsing, d::Diagnostics)
    lsize = link_list_update!(rng, qmc_state, H, d)
    return multibranch_cluster_update!(rng, lsize, qmc_state, H, d)
end
