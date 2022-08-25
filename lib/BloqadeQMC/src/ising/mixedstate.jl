
########################## finite-beta #######################################

function mc_step_beta!(f::Function, rng::AbstractRNG, qmc_state::BinaryThermalState, H::AbstractIsing, beta::Real, d::Diagnostics; eq::Bool = false, kw...)
    num_ops = full_diagonal_update_beta!(rng, qmc_state, H, beta; eq=eq)
    lsize = cluster_update!(rng, qmc_state, H, d; kw...)
    f(lsize, qmc_state, H)
    return num_ops
end

mc_step_beta!(f::Function, qmc_state, H, beta, d::Diagnostics; eq = false, kw...) = mc_step_beta!(f, Random.GLOBAL_RNG, qmc_state, H, beta, d; eq = eq, kw...)
mc_step_beta!(rng::AbstractRNG, qmc_state, H, beta, d::Diagnostics; eq = false, kw...) = mc_step_beta!((args...) -> nothing, rng, qmc_state, H, beta, d; eq = eq, kw...)
mc_step_beta!(qmc_state, H, beta, d::Diagnostics; eq = false, kw...) = mc_step_beta!(Random.GLOBAL_RNG, qmc_state, H, beta, d; eq = eq, kw...)

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


function full_diagonal_update_beta!(rng::AbstractRNG, qmc_state::BinaryThermalState, H::AbstractIsing, beta::Real; eq::Bool = false)
    P_norm = beta * diag_update_normalization(H)

    num_ids = count(isidentity(H), qmc_state.operator_list)

    spin_prop = copyto!(qmc_state.propagated_config, qmc_state.left_config)  # the propagated spin state

    @inbounds for (n, op) in enumerate(qmc_state.operator_list)
        if !isdiagonal(H, op)
            spin_prop[getsite(H, op)] ⊻= 1  # spinflip
        elseif !isidentity(H, op)
            if rand(rng)*P_norm < (num_ids + 1)
                qmc_state.operator_list[n] = makeidentity(H)
                num_ids += 1
            end
        else
            if rand(rng)*num_ids < P_norm
                op, _ = insert_diagonal_operator!(rng, qmc_state, H, spin_prop, n)
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
        resize_op_list!(qmc_state, H, round(Int, 1.5*num_ops, RoundUp))
    end

    return num_ops
end
full_diagonal_update_beta!(qmc_state, H, beta; eq = false) = full_diagonal_update_beta!(Random.GLOBAL_RNG, qmc_state, H, beta; eq = eq)
