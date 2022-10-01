# Todo: Merge the following two mc_step functions

# mc_step from groundstate.jl

function mc_step!(f::Function, rng::AbstractRNG, qmc_state::BinaryGroundState, H::Hamiltonian, d::Diagnostics; kw...)
    full_diagonal_update!(rng, qmc_state, H, d)
    lsize = cluster_update!(rng, qmc_state, H, d; kw...)
    f(lsize, qmc_state, H)
end
# mc_step!(f::Function, qmc_state::BinaryGroundState, H::Hamiltonian, d::Diagnostics; kw...) = mc_step!(f, Random.GLOBAL_RNG, qmc_state, H, d; kw...)
mc_step!(rng::AbstractRNG, qmc_state::BinaryGroundState, H::Hamiltonian, d::Diagnostics; kw...) = mc_step!((args...) -> nothing, rng, qmc_state, H, d; kw...)
# mc_step!(qmc_state::BinaryGroundState, H::Hamiltonian, d::Diagnostics; kw...) = mc_step!(Random.GLOBAL_RNG, qmc_state, H, d; kw...)


# mc_step from mixedstate.jl

function mc_step_beta!(f::Function, rng::AbstractRNG, qmc_state::BinaryThermalState, H::AbstractIsing, beta::Real, d::Diagnostics; eq::Bool = false, kw...)
    num_ops = full_diagonal_update_beta!(rng, qmc_state, H, beta; eq=eq)
    lsize = cluster_update!(rng, qmc_state, H, d; kw...)
    f(lsize, qmc_state, H)
    return num_ops
end

# mc_step_beta!(f::Function, qmc_state, H, beta, d::Diagnostics; eq = false, kw...) = mc_step_beta!(f, Random.GLOBAL_RNG, qmc_state, H, beta, d; eq = eq, kw...)
mc_step_beta!(rng::AbstractRNG, qmc_state, H, beta, d::Diagnostics; eq = false, kw...) = mc_step_beta!((args...) -> nothing, rng, qmc_state, H, beta, d; eq = eq, kw...)
# mc_step_beta!(qmc_state, H, beta, d::Diagnostics; eq = false, kw...) = mc_step_beta!(Random.GLOBAL_RNG, qmc_state, H, beta, d; eq = eq, kw...)