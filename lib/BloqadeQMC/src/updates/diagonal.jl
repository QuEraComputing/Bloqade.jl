# deleted alignment_check and insert_diagonal_operator for TFIM case

function insert_diagonal_operator!(rng::AbstractRNG, qmc_state::BinaryQMCState{K, V}, H::AbstractLTFIM{<:AbstractImprovedOperatorSampler{K, T}}, spin_prop::V, n::Int) where {K, T, V}
# Gets called in full_diagonal_update() below
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

# Diagonal update from groundstate.jl

function full_diagonal_update!(rng::AbstractRNG, qmc_state::BinaryGroundState, H::AbstractIsing, d::Diagnostics)
# gets called in mc_step() in mc_step.jl
    spin_prop = copyto!(qmc_state.propagated_config, qmc_state.left_config)  # the propagated spin state
    runstats = d.runstats

    if !(runstats isa NoStats) # Technically, this runstats stuff could be cut, but it is useful for checking MC on developer side
        failures = 0
        count = 0
    end

    for (n, op) in enumerate(qmc_state.operator_list)
        if !isdiagonal(H, op)
            @inbounds spin_prop[getsite(H, op)] ⊻= 1  # spinflip
        else
            op = nothing
            if !(runstats isa NoStats); i = -1; end
            while op === nothing # here we are attempting step 3 until we achieve a successful insertion
                op, _ = insert_diagonal_operator!(rng, qmc_state, H, spin_prop, n)
                if !(runstats isa NoStats); i += 1; end
            end
            if !(runstats isa NoStats); failures += i; count += 1; end
        end
    end

    # @debug statements are not run unless debug logging is enabled
    @debug("Diagonal Update basis state propagation status: $(spin_prop == qmc_state.right_config)",
           spin_prop,
           qmc_state.right_config)

    if !(runstats isa NoStats)
        fit!(runstats, :diag_update_fails, failures/count)
    end
end

#########################################################################################################################################

# Diagonal updates (full_diagonal_update_beta and resize_op_list from mixedstate.jl)

function resize_op_list!(qmc_state::BinaryThermalState{K}, H::AbstractIsing, new_size::Int) where {K}
    operator_list = filter!(!isidentity(H), qmc_state.operator_list)
    len = length(operator_list)

    if len < new_size
        tail = init_op_list(new_size - len, K)
        append!(operator_list, tail)
    end

    len = 4*length(operator_list)
    # these are going to be overwritten by the cluster update which will be
    # called right after the diagonal update that called this function
    resize!(qmc_state.linked_list, len)
    resize!(qmc_state.leg_types, len)
    resize!(qmc_state.associates, len)
    resize!(qmc_state.leg_sites, len)
    resize!(qmc_state.op_indices, len)
    resize!(qmc_state.in_cluster, len)
end


function full_diagonal_update_beta!(rng::AbstractRNG, qmc_state::BinaryThermalState, H::AbstractIsing, beta::Real; eq::Bool = false) # eq = true must be set during eq phase for tuning M (resizing)
# gets called in mc_step() in mc_step.jl
    # As we can see: thermal case does not call runstats
    P_norm = beta * diag_update_normalization(H)

    num_ids = count(isidentity(H), qmc_state.operator_list)

    spin_prop = copyto!(qmc_state.propagated_config, qmc_state.left_config)  # the propagated spin state

    @inbounds for (n, op) in enumerate(qmc_state.operator_list)
        if !isdiagonal(H, op) # Step 4
            spin_prop[getsite(H, op)] ⊻= 1  # spinflip
        elseif !isidentity(H, op)
            if rand(rng)*P_norm < (num_ids + 1) # implements condition (16)  
                qmc_state.operator_list[n] = makeidentity(H) # This corresponds to removing the operator in Step 1
                num_ids += 1
            end
        else
            if rand(rng)*num_ids < P_norm # implements condition (17)
                op, _ = insert_diagonal_operator!(rng, qmc_state, H, spin_prop, n) # This corresponds to inserting an operator in Step 2. 
                # Sampling step 3 is implemented in insert_diagonal_operator!()
                if op !== nothing
                    num_ids -= 1
                end
            end
        end
    end

    # DEBUG
    # if spin_prop != qmc_state.right_config  # check the spin propagation for error
    #     error("Basis state propagation error in diagonal update!")
    # end

    total_list_size = length(qmc_state.operator_list)
    num_ops = total_list_size - num_ids

    if eq && 1.2*num_ops > total_list_size
        resize_op_list!(qmc_state, H, round(Int, 1.5*num_ops, RoundUp)) # Only happens in equilibration phase
    end

    return num_ops
end

