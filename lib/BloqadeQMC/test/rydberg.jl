using Test, BloqadeQMC
using Statistics
using Random
using RandomNumbers
using Measurements
using Measurements: value, uncertainty
using BinningAnalysis
using BloqadeQMC: Chain, Square
using BloqadeLattices: generate_sites, ChainLattice
using BloqadeExpr: rydberg_h 
using Yao: mat, ArrayReg
using LinearAlgebra
#using PLots

# Generate ED values - do we want the ED to run every time? Or do we want to pre-calculate and store the values in a dict?

βs = [0.005, 0.05, 0.5]
nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.48)

Ω = 2π * 4
Δ_step = 30
Δ = LinRange(-2π * 9, 2π * 9, Δ_step)

energy_ED = zeros(3, Δ_step)

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
    h_m = Matrix(mat(h_ii)) 
    energies, vecs = LinearAlgebra.eigen(h_m) 
    N_e = length(energies)

    for (jj, β) in enumerate(βs)
        w = exp.(-β .* (energies .- energies[1]))
        energy_ED[jj, ii] = sum(w .* energies) / sum(w)
    end
end

### Now, start running QMC tests.

# THRESHOLD_t = 2.576             threshold for t-test with ∞ DOF and 99.5% confidence   
THRESHOLD_χ = 43.77             # threshold for χ² test with 30 DOF and p=0.05

R_b = 7.74                      # This does correspond to C6 = 862690 * 2pi Mhz μm^6 at a Rabi frequency of Ω = 2π * 4.
Ω = 2π * 4
N = 9
a = 5.48

lat = Chain(N, a, false; trunc=Inf)

@testset "1D Chain (9 atoms), β=0.005" begin
    β = 0.005
    EQ_MCS = 30
    MCS = 100_000
    M = 50

    rng = MersenneTwister(3214)

    energy_QMC_β1 = []

    χ_squared = 0

    for ii in 1:Δ_step
        H = Rydberg(lat, R_b, Ω, Δ[ii])
        ts = BinaryThermalState(H, M)
        d = Diagnostics()
    
        [mc_step_beta!(rng, ts, H, β, d, eq=true) for i in 1:EQ_MCS] #equilibration phase
    
        ns = zeros(MCS)
    
        for i in 1:MCS # Monte Carlo Steps
            ns[i] = mc_step_beta!(rng, ts, H, β, d, eq=false)
        end

        append!(energy_QMC_β1, mean_and_stderr(x -> -x/β, ns) + H.energy_shift)
        χ_squared += abs2(value(energy_QMC_β1[ii]) - energy_ED[1, ii]) / abs2(uncertainty(energy_QMC_β1[ii]))

        # @test abs(stdscore(energy_QMC_β1[ii], energy_ED[1,ii])) < THRESHOLD_t       t-test for testing each QMC run individually
    end
    @test χ_squared < THRESHOLD_χ
end


@testset "1D Chain (9 atoms), β=0.05" begin
    β = 0.05
    EQ_MCS = 100
    MCS = 100_000
    M = 500

    rng = MersenneTwister(13542)

    energy_QMC_β2 = []

    χ_squared = 0

    for ii in 1:Δ_step
        H = Rydberg(lat, R_b, Ω, Δ[ii])
        ts = BinaryThermalState(H, M)
        d = Diagnostics()
    
        [mc_step_beta!(rng, ts, H, β, d, eq=true) for i in 1:EQ_MCS] #equilibration phase
    
        ns = zeros(MCS)
    
        for i in 1:MCS # Monte Carlo Steps
            ns[i] = mc_step_beta!(rng, ts, H, β, d, eq=false)
        end

        append!(energy_QMC_β2, mean_and_stderr(x -> -x/β, ns) + H.energy_shift)
        χ_squared += abs2(value(energy_QMC_β2[ii]) - energy_ED[2, ii]) / abs2(uncertainty(energy_QMC_β2[ii]))
        # @test abs(stdscore(energy_QMC_β2[ii], energy_ED[2,ii])) < THRESHOLD_t
    end
    @test χ_squared < THRESHOLD_χ
end

@testset "1D Chain (9 atoms), β=0.5" begin
    β = 0.5
    EQ_MCS = 100
    MCS = 10_000
    M = 5000

    rng = MersenneTwister(12345)

    energy_QMC_β3 = []

    χ_squared = 0

    for ii in 1:Δ_step
        H = Rydberg(lat, R_b, Ω, Δ[ii])
        ts = BinaryThermalState(H, M)
        d = Diagnostics()
    
        [mc_step_beta!(rng, ts, H, β, d, eq=true) for i in 1:EQ_MCS] #equilibration phase
        
        ns = zeros(MCS)
    
        for i in 1:MCS # Monte Carlo Steps
            ns[i] = mc_step_beta!(rng, ts, H, β, d, eq=false)
        end

        append!(energy_QMC_β3, mean_and_stderr(x -> -x/β, ns) + H.energy_shift)

        χ_squared += abs2(value(energy_QMC_β3[ii]) - energy_ED[3, ii]) / abs2(uncertainty(energy_QMC_β3[ii]))

        # @test abs(stdscore(energy_QMC_β3[ii], energy_ED[3,ii])) < THRESHOLD_t
    end
    @test χ_squared < THRESHOLD_χ
end