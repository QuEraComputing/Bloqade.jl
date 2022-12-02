module BloqadeQMC

using Measurements
using Statistics
import Statistics: var, mean
using Distributions
using FFTW

using DelimitedFiles
using JLD2
using Printf

using PushVectors
import PushVectors: PushVector
using DataStructures
using SparseArrays
using LinearAlgebra

using Random

import BloqadeLattices: rydberg_interaction_matrix, AtomList
import BloqadeExpr: RydbergHamiltonian, get_rydberg_params, is_time_function

import Base: zero, one, convert
import Base: length, size, eltype, setindex!, getindex, firstindex, lastindex
import Base: rand, show, pop!, push!, append!, isempty, empty!, count

using BinningAnalysis
import BinningAnalysis: varN, std_error

export BinaryQMCState, BinaryGroundState, BinaryThermalState,
        Hamiltonian, AbstractIsing, AbstractTFIM, TFIM, AbstractLTFIM, LTFIM, GeneralLTFIM, AbstractRydberg, Rydberg,
        haslongitudinalfield,
        nspins, nbonds, energy, energy_density, mc_step!, mc_step_beta!,
        resize_op_list!,
        sample, simulation_cell, magnetization, staggered_magnetization, domain_wall_density, kagome_nematic, correlation_functions,
        num_single_site_diag, num_single_site_offdiag, update_two_pt_fn!,
        num_single_site, num_two_site_diag, autocorrelation, correlation_time, jackknife, bootstrap, mean_and_stderr,
        Lattice, OneDLattice, BravaisLattice, PolyLattice,
        Triangle, Rectangle, Kagome, Ruby, Custom, lattice_bond_spins, distance_matrix,
        ProbabilityAlias, probability_vector

export AbstractRunStats, NoStats, RunStats, RunStatsHistogram, Diagnostics
export NoTransitionMatrix, TransitionMatrix, PositionDependentTransitionMatrix, CombinedTransitionMatrix
export ProductState, PlusState, AbstractProductState, AbstractTrialState

export rydberg_qmc

@inline function pop!(v::PushVector)
    @boundscheck isempty(v) && throw(ArgumentError("vector must be non-empty"))
    x = @inbounds v.parent[v.len]
    v.len -= 1
    x
end


include("lattices/bond_spins.jl")
include("lattices/Lattices.jl")
include("probabilityvectors/probabilityvector.jl")
include("operatorsamplers/operatorsampler.jl")
include("trialstates/trialstate.jl")
include("qmc_state.jl")
include("hamiltonian.jl")
include("operatorsamplers/improved_op_sampler.jl")
include("diagnostics/Diagnostics.jl")
include("ising/Ising.jl")
include("measurements.jl")
include("error.jl")


end
