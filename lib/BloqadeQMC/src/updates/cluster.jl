# Main parts here: huge link_list_update!() and cluster_update!() functions.
# They don´t seem to call each other -> understand when each is used & how they interact.

function link_list_update!(::AbstractRNG, qmc_state::BinaryQMCState, H::AbstractIsing, ::Diagnostics) # fairly standard, builds necessary data structure for both line and multibranch 
    Ns = nspins(H)
    spin_left = qmc_state.left_config

    # retrieve linked list data structures
    LinkList = qmc_state.linked_list  # needed for cluster update
    LegType = qmc_state.leg_types

    # A diagonal bond operator has non trivial associates for cluster building
    Associates = qmc_state.associates

    op_indices = qmc_state.op_indices
    leg_sites = qmc_state.leg_sites

    if qmc_state isa BinaryGroundState
        First = qmc_state.first

        # The first N elements of the linked list are the spins of the LHS basis state
        @inbounds for i in 1:Ns
            LegType[i] = spin_left[i]
            Associates[i] = 0
            First[i] = i
            if H isa AbstractLTFIM
                op_indices[i] = 0
                leg_sites[i] = 0
            end
        end
        idx = Ns
    else
        First = fill!(qmc_state.first, 0)  #initialize the First list
        Last = fill!(qmc_state.last, 0)   #initialize the Last list
        idx = 0
    end

    spin_prop = copyto!(qmc_state.propagated_config, spin_left)  # the propagated spin state

    # Now, add the 2M operators to the linked list. Each has either 2 or 4 legs
    @inbounds for (n, op) in enumerate(qmc_state.operator_list)
        if issiteoperator(H, op)
            site = getsite(H, op)
            # lower or left leg
            idx += 1
            F = First[site]
            LinkList[idx] = F
            if qmc_state isa BinaryGroundState || !iszero(F)
                LinkList[F] = idx  # completes backwards link
            else
                Last[site] = idx
            end

            LegType[idx] = spin_prop[site]
            Associates[idx] = 0
            if H isa AbstractLTFIM
                op_indices[idx] = n
                leg_sites[idx] = 1
            end

            if !isdiagonal(H, op)  # off-diagonal site operator
                spin_prop[site] ⊻= 1  # spinflip
            end

            # upper or right leg
            idx += 1
            First[site] = idx
            LegType[idx] = spin_prop[site]
            Associates[idx] = 0
            if H isa AbstractLTFIM
                op_indices[idx] = n
                leg_sites[idx] = 1
            end
        elseif qmc_state isa BinaryGroundState || isbondoperator(H, op)  # diagonal bond operator
            site1, site2 = bond = getbondsites(H, op)
            spins = spin_prop[site1], spin_prop[site2]
            num_sites = 2  # length(bond)
            num_legs = 4   # 2*num_sites

            # vertex leg indices
            # thinking of imaginary time as increasing as we move upward,
            # these indices refer to the
            # lower left, lower right, upper left, upper right legs respectively
            # v1, v2, v3, v4 = idx + 1, idx + 2, idx + 3, idx + 4

            @simd for i in 1:num_sites
                v = idx + i
                st = bond[i]
                F = First[st]
                LinkList[v] = F
                if qmc_state isa BinaryGroundState || !iszero(F)
                    LinkList[F] = v  # completes backwards link
                else
                    Last[st] = v
                end
                First[st] = v + num_sites
            end

            @simd for i in 1:num_legs
                v = idx + i
                m = mod(i, 1:num_sites)
                LegType[v] = spins[m]
                Associates[v] = v + 1
                if H isa AbstractLTFIM
                    op_indices[v] = n
                    leg_sites[v] = m
                end
            end
            Associates[idx + num_legs] = idx + 1
            idx += num_legs
        end
    end

    if qmc_state isa BinaryGroundState
        # The last N elements of the linked list are the final spin state
        @inbounds for i in 1:Ns
            idx += 1
            F = First[i]
            LinkList[idx] = F
            LinkList[F] = idx
            LegType[idx] = spin_prop[i]
            Associates[idx] = 0
            if H isa AbstractLTFIM
                op_indices[idx] = 0
                leg_sites[idx] = 0
            end
        end
    else
        #Periodic boundary conditions for finite-beta
        @inbounds for i in 1:Ns
            F = First[i]
            if !iszero(F)  #This might be encountered at high temperatures
                L = Last[i]
                LinkList[F] = L
                LinkList[L] = F
            end
        end
    end
    # @debug statements are not run unless debug logging is enabled
    @debug("Link List basis state propagation status: $(spin_prop == qmc_state.right_config)",
           spin_prop,
           qmc_state.right_config)

    return idx
end
link_list_update!(qmc_state, H, d::Diagnostics) =
    link_list_update!(Random.GLOBAL_RNG, qmc_state, H, d)


#############################################################################


@inline function _map_back_operator_list!(ocount::Int, qmc_state::BinaryQMCState, H::AbstractIsing, d::Diagnostics)
    operator_list = qmc_state.operator_list
    LegType = qmc_state.leg_types

    # if we build an array that maps n to ocount, this loop will become
    #   very easy to parallelize
    @inbounds for (n, op) in enumerate(operator_list)
        if isbondoperator(H, op)
            if H isa AbstractLTFIM
                s1, s2 = LegType[ocount], LegType[ocount+1]
                t = getbondtype(H, s1, s2)
                operator_list[n] = convertoperatortype(H, op, t)
                fit!(d.tmatrix, op, operator_list[n])
            end
            ocount += 4
        elseif issiteoperator(H, op)
            if LegType[ocount] == LegType[ocount+1]  # diagonal
                operator_list[n] = makediagonalsiteop(H, getsite(H, op))
            else  # off-diagonal
                operator_list[n] = makeoffdiagonalsiteop(H, getsite(H, op))
            end
            fit!(d.tmatrix, op, operator_list[n])
            ocount += 2
        end
    end
end


@inline function _map_back_basis_states!(rng::AbstractRNG, lsize::Int, qmc_state::BinaryGroundState, H::AbstractIsing)
    Ns = nspins(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config
    LegType = qmc_state.leg_types

    @inbounds for i in 1:Ns
        spin_left[i] = LegType[i]  # left basis state
        spin_right[i] = LegType[lsize-Ns+i]  # right basis state
    end
    return Ns + 1  # next one is leg Ns + 1
end

@inline function _map_back_basis_states!(rng::AbstractRNG, lsize::Int, qmc_state::BinaryThermalState, H::AbstractIsing)
    Ns = nspins(H)
    spin_left, spin_right = qmc_state.left_config, qmc_state.right_config
    LegType = qmc_state.leg_types

    First = qmc_state.first
    Last = qmc_state.last
    @inbounds for i in 1:Ns
        F = First[i]
        if !iszero(F)
            spin_left[i] = LegType[Last[i]]  # left basis state
            spin_right[i] = LegType[F]  # right basis state
        else
            #randomly flip spins not connected to operators
            spin_left[i] = spin_right[i] = rand(rng, Bool)
        end
    end
    return 1  # first leg
end


function trialstate_weight_change(qmc_state::BinaryGroundState, lsize::Int, Ns::Int, i::Int) # This we can get rid of because it's based on TrialState, only returns logweights -> before cutting this, just let it return 0
    if !(qmc_state.trialstate isa AbstractProductState)
        if i <= Ns
            push!(qmc_state.trialstate.left_flips, i)
        elseif i > (lsize - Ns)
            push!(qmc_state.trialstate.right_flips, i - lsize + Ns)
        end
        return 0.0
    else
        if i <= Ns
            # return logweightchange(qmc_state.trialstate, qmc_state.left_config[i])
            return 0.0
        elseif i > (lsize - Ns)
            # return logweightchange(qmc_state.trialstate, qmc_state.right_config[i])
            return 0.0
        else
            return 0.0
        end
    end
end
trialstate_weight_change(qmc_state::BinaryThermalState, lsize::Int, Ns::Int, i::Int) = 0.0

#############################################################################

# changed the function name below from cluster_update!() to generic_cluster_update!() in order to distinguish it from the other cluster_update function defined at the very bottom

function generic_cluster_update!(rng::AbstractRNG, update_kernel!::Function, acceptance::Function, lsize::Int, qmc_state::BinaryQMCState, H::AbstractIsing, d::Diagnostics)
    Ns = nspins(H)
    operator_list = qmc_state.operator_list

    LinkList = qmc_state.linked_list
    LegType = qmc_state.leg_types
    Associates = qmc_state.associates
    op_indices = qmc_state.op_indices

    in_cluster = fill!(qmc_state.in_cluster, 0)
    cstack = qmc_state.cstack # This is the stack of vertices in a cluster
    current_cluster = qmc_state.current_cluster
    runstats = d.runstats

    ccount = 0  # cluster number counter

    num_accept = 0
    num_reject = 0

    @inbounds for i in 1:lsize
        # Add a new leg onto the cluster
        if iszero(in_cluster[i]) && iszero(Associates[i])
            ccount += 1
            push!(cstack, i)
            in_cluster[i] = ccount
            if qmc_state isa BinaryGroundState && !(qmc_state.trialstate isa AbstractProductState)
                left_flips = empty!(qmc_state.trialstate.left_flips)
                right_flips = empty!(qmc_state.trialstate.right_flips)
            end

            empty!(current_cluster)
            push!(current_cluster, i)
            lnA = 0.0
            if qmc_state isa BinaryGroundState
                lnA += trialstate_weight_change(qmc_state, lsize, Ns, i)
            end

            while !isempty(cstack)
                leg = LinkList[pop!(cstack)]

                if iszero(in_cluster[leg])
                    in_cluster[leg] = ccount  # add the new leg and flip it
                    push!(current_cluster, leg)
                    if qmc_state isa BinaryGroundState
                        lnA += trialstate_weight_change(qmc_state, lsize, Ns, i)
                    end
                    a = Associates[leg]

                    a == 0 && continue
                    # from this point on, we know we're on a bond op
                    op = operator_list[op_indices[leg]]
                    w = getweightindex(H, op) - getoperatortype(H, op)
                    preflip_bond_type, postflip_bond_type = update_kernel!(qmc_state, H, ccount, leg, a)
                    lnA += (
                        getlogweight(H.op_sampler, w + postflip_bond_type)
                        - getlogweight(H.op_sampler, w + preflip_bond_type)
                    )
                end
            end

            if qmc_state isa BinaryGroundState && !(qmc_state.trialstate isa AbstractProductState)
                lnA += logweightchange(qmc_state.trialstate, qmc_state.left_config, left_flips)
                lnA += logweightchange(qmc_state.trialstate, qmc_state.right_config, right_flips)
            end

            A = acceptance(H, lnA)
            fit!(runstats, :cluster_update_accept, A)

            if rand(rng) < A
                @inbounds for j in current_cluster
                    LegType[j] ⊻= 1  # spinflip
                end
                fit!(runstats, :accepted_cluster_sizes, length(current_cluster))
                num_accept += 1
            else
                fit!(runstats, :rejected_cluster_sizes, length(current_cluster))
                num_reject += 1
            end
            fit!(runstats, :cluster_sizes, length(current_cluster))
        end
    end

    if !(runstats isa NoStats)
        fit!(runstats, :accepted_cluster_count, num_accept+1)
        fit!(runstats, :rejected_cluster_count, num_reject+1)
        fit!(runstats, :cluster_count, ccount+1)
    end

    # map back basis states and operator list
    ocount = _map_back_basis_states!(rng, lsize, qmc_state, H)
    _map_back_operator_list!(ocount, qmc_state, H, d)

    return lsize
end

# might put this cluster_update function in a more visible position somewhere

function cluster_update!(rng, qmc_state, H::Hamiltonian, d::Diagnostics; p::Float64=0.0, kw...)
    if rand(rng) < p
        multibranch_update!(rng, qmc_state, H, d)
    else
        line_update!(rng, qmc_state, H, d)
    end
end
