# cluster_update!(rng, qmc_state, H::Hamiltonian, runstats; kw...) = multibranch_update!(rng, qmc_state, H, runstats)
function cluster_update!(rng, qmc_state, H::Hamiltonian, d::Diagnostics; p::Float64=0.0, kw...)
    if rand(rng) < p
        multibranch_update!(rng, qmc_state, H, d)
    else
        line_update!(rng, qmc_state, H, d)
    end
end

function mc_step!(f::Function, rng::AbstractRNG, qmc_state::BinaryGroundState, H::Hamiltonian, d::Diagnostics; kw...)
    full_diagonal_update!(rng, qmc_state, H, d)
    lsize = cluster_update!(rng, qmc_state, H, d; kw...)
    f(lsize, qmc_state, H)
end
mc_step!(f::Function, qmc_state::BinaryGroundState, H::Hamiltonian, d::Diagnostics; kw...) = mc_step!(f, Random.GLOBAL_RNG, qmc_state, H, d; kw...)
mc_step!(rng::AbstractRNG, qmc_state::BinaryGroundState, H::Hamiltonian, d::Diagnostics; kw...) = mc_step!((args...) -> nothing, rng, qmc_state, H, d; kw...)
mc_step!(qmc_state::BinaryGroundState, H::Hamiltonian, d::Diagnostics; kw...) = mc_step!(Random.GLOBAL_RNG, qmc_state, H, d; kw...)


function full_diagonal_update!(rng::AbstractRNG, qmc_state::BinaryGroundState, H::AbstractIsing, d::Diagnostics)
    spin_prop = copyto!(qmc_state.propagated_config, qmc_state.left_config)  # the propagated spin state
    runstats = d.runstats

    if !(runstats isa NoStats)
        failures = 0
        count = 0
    end

    for (n, op) in enumerate(qmc_state.operator_list)
        if !isdiagonal(H, op)
            @inbounds spin_prop[getsite(H, op)] ⊻= 1  # spinflip
        else
            op = nothing
            if !(runstats isa NoStats); i = -1; end
            while op === nothing
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
full_diagonal_update!(qmc_state, H, d::Diagnostics) = full_diagonal_update!(Random.GLOBAL_RNG, qmc_state, H, runstats)
