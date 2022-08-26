# TODO: there might be a slight performance degradation with this
#       functional approach, might wanna try to inject this function directly
#       into the body of the cluster update
@inline function line_kernel!(qmc_state::BinaryQMCState, H::AbstractIsing, ccount::Int, leg::Int, a::Int)
    Ns = nspins(H)
    LegType, Associates, leg_sites = qmc_state.leg_types, qmc_state.associates, qmc_state.leg_sites
    in_cluster, cstack, current_cluster = qmc_state.in_cluster, qmc_state.cstack, qmc_state.current_cluster

    @inbounds ll, la = LegType[leg], LegType[a]
    @inbounds sl, sa = leg_sites[leg], leg_sites[a]
    # TODO: check if this is inputting the spins in the correct order
    if sl > sa
        preflip_bond_type = getbondtype(H, ll, la)
        postflip_bond_type = getbondtype(H, !ll, la)
    else
        preflip_bond_type = getbondtype(H, la, ll)
        postflip_bond_type = getbondtype(H, la, !ll)
    end
    # ^using circshift for this is too slow, gonna have to manually expand out
    #  a bunch of ifs for k-local (might need to do some metaprogramming!)
    # using SVectors and sortperm should work (it's super fast!) and is easier
    # also, should short-circuit if order doesn't matter, i.e.
    #  the longitudinal field is uniform + coordination numbers are all equal

    # now add the straight-through leg to the cluster
    @inbounds straight_thru = Associates[a]
    @inbounds in_cluster[straight_thru] = ccount
    push!(cstack, straight_thru)
    push!(current_cluster, straight_thru)

    return preflip_bond_type, postflip_bond_type
end

line_acceptance(H::AbstractIsing, lnA::T) where {T <: Real} = exp(min(lnA, zero(T)))

function line_cluster_update!(rng::AbstractRNG, lsize::Int, qmc_state::BinaryQMCState, H::AbstractIsing, d::Diagnostics)
    return cluster_update!(rng, line_kernel!, line_acceptance, lsize, qmc_state, H, d)
end

line_cluster_update!(lsize, qmc_state, H, d::Diagnostics) = line_cluster_update!(Random.GLOBAL_RNG, lsize, qmc_state, H, d)


#############################################################################


function line_update!(rng::AbstractRNG, qmc_state::BinaryQMCState, H::AbstractIsing, d::Diagnostics)
    lsize = link_list_update!(rng, qmc_state, H, d)
    return line_cluster_update!(rng, lsize, qmc_state, H, d)
end
line_update!(qmc_state::BinaryQMCState, H::AbstractIsing, d::Diagnostics) = line_update!(Random.GLOBAL_RNG, qmc_state, H, d)
