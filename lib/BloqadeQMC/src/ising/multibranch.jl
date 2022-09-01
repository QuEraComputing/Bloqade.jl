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
    return cluster_update!(rng, multibranch_kernel!, multibranch_acceptance, lsize, qmc_state, H, d)
end

# simplified implementation for TFIM
function multibranch_cluster_update!(rng::AbstractRNG, lsize::Int, qmc_state::BinaryQMCState, H::AbstractTFIM, d::Diagnostics)
    Ns = nspins(H)
    operator_list = qmc_state.operator_list

    LinkList = qmc_state.linked_list
    LegType = qmc_state.leg_types
    Associates = qmc_state.associates

    in_cluster = fill!(qmc_state.in_cluster, false)
    cstack = qmc_state.cstack  # This is the stack of vertices in a cluster
    runstats = d.runstats

    if !(runstats isa NoStats)
        ccount = 0  # cluster number counter
    end

    @inbounds for i in 1:lsize
        # Add a new leg onto the cluster
        if (iszero(in_cluster[i]) && Associates[i] == 0)
            if !(runstats isa NoStats)
                cluster_size = 1
                ccount += 1
            end
            push!(cstack, i)
            in_cluster[i] = true

            flip = rand(rng, Bool)  # flip a coin for the SW cluster flip
            if flip
                LegType[i] ⊻= 1  # spinflip
            end

            while !isempty(cstack)
                leg = LinkList[pop!(cstack)]

                if iszero(in_cluster[leg])
                    in_cluster[leg] = true  # add the new leg and flip it
                    if !(runstats isa NoStats); cluster_size += 1; end
                    if flip
                        LegType[leg] ⊻= 1
                    end

                    # now check all associates and add them to the cluster
                    a = Associates[leg]
                    while a != 0 && iszero(in_cluster[a])
                        push!(cstack, a)
                        in_cluster[a] = true
                        if !(runstats isa NoStats); cluster_size += 1; end
                        if flip
                            LegType[a] ⊻= 1
                        end
                        a = Associates[a]
                    end
                end
            end
            if !(runstats isa NoStats); fit!(runstats, :cluster_sizes, float(cluster_size)); end
            fit!(runstats, :cluster_update_accept, 1/2)
        end
    end

    if !(runstats isa NoStats)
        fit!(runstats, :cluster_count, float(ccount))
    end

    # map back basis states and operator list
    ocount = _map_back_basis_states!(rng, lsize, qmc_state, H)
    _map_back_operator_list!(ocount, qmc_state, H, d)

    return lsize
end
multibranch_cluster_update!(lsize, qmc_state, H, d::Diagnostics) =
    multibranch_cluster_update!(Random.GLOBAL_RNG, lsize, qmc_state, H, d)


#############################################################################


function multibranch_update!(rng::AbstractRNG, qmc_state::BinaryQMCState, H::AbstractIsing, d::Diagnostics)
    lsize = link_list_update!(rng, qmc_state, H, d)
    return multibranch_cluster_update!(rng, lsize, qmc_state, H, d)
end
multibranch_update!(qmc_state::BinaryQMCState, H::AbstractIsing, d::Diagnostics) =
    multibranch_update!(Random.GLOBAL_RNG, qmc_state, H, d)
