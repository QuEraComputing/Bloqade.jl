using Test
using Statistics
using Random
using RandomNumbers
using Measurements
using Measurements: value, uncertainty
using BloqadeLattices: generate_sites, ChainLattice
using BloqadeQMC: rydberg_QMC, BinaryThermalState, Diagnostics, mc_step_beta!
# using BloqadeQMC
using BloqadeExpr: rydberg_h 
using Yao: mat, ArrayReg
using LinearAlgebra
using BinningAnalysis

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

@testset "1D Chain (9 atoms), β=0.005" begin
    β = 0.005
    EQ_MCS = 30
    MCS = 100_000
    M = 50

    rng = MersenneTwister(3214)

    energy_QMC_β1 = []

    χ_squared = 0

    for ii in 1:Δ_step
        @show ii
        h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
        H = rydberg_QMC(h_ii)
        ts = BinaryThermalState(H, M)
        d = Diagnostics()
    
        [mc_step_beta!(rng, ts, H, β, d, eq=true) for i in 1:EQ_MCS] #equilibration phase
    
        ns = zeros(MCS)
    
        for i in 1:MCS # Monte Carlo Steps
            ns[i] = mc_step_beta!(rng, ts, H, β, d, eq=false)
        end

        # Binning analysis 
        energy(x) = -x / β + H.energy_shift
        BE = LogBinner(energy.(ns))
        τ_energy = tau(BE)
        ratio = 2 * τ_energy + 1
        energy_binned = measurement(mean(BE), std_error(BE)*sqrt(ratio)) 
        append!(energy_QMC_β1, energy_binned)

        χ_squared += abs2(value(energy_QMC_β1[ii]) - energy_ED[1, ii]) / abs2(uncertainty(energy_QMC_β1[ii]))
        @show χ_squared
        # @test abs(stdscore(energy_QMC_β1[ii], energy_ED[1,ii])) < THRESHOLD_t       t-test for testing each QMC run individually
    end

    @test χ_squared < THRESHOLD_χ

    # scatter(Δ/2π, value.(energy_QMC_β1); yerror=uncertainty.(energy_QMC_β1), marker=:x)
    # scatter!(Δ/2π, energy_ED[1,:], marker=:x)
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
        @show ii
        h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
        H = rydberg_QMC(h_ii)
        ts = BinaryThermalState(H, M)
        d = Diagnostics()
    
        [mc_step_beta!(rng, ts, H, β, d, eq=true) for i in 1:EQ_MCS] #equilibration phase
    
        ns = zeros(MCS)
    
        for i in 1:MCS # Monte Carlo Steps
            ns[i] = mc_step_beta!(rng, ts, H, β, d, eq=false)
        end

        # Binning analysis 
        energy(x) = -x / β + H.energy_shift
        BE = LogBinner(energy.(ns))
        τ_energy = tau(BE)
        ratio = 2 * τ_energy + 1
        energy_binned = measurement(mean(BE), std_error(BE)*sqrt(ratio)) 
        append!(energy_QMC_β2, energy_binned)
        println()
        #append!(energy_QMC_β2, mean_and_stderr(x -> -x/β, ns) + H.energy_shift)

        @show χ_squared += abs2(value(energy_QMC_β2[ii]) - energy_ED[2, ii]) / abs2(uncertainty(energy_QMC_β2[ii]))
        # @test abs(stdscore(energy_QMC_β2[ii], energy_ED[2,ii])) < THRESHOLD_t
    end
    @test χ_squared < THRESHOLD_χ

    # scatter(Δ/2π, value.(energy_QMC_β2); yerror=uncertainty.(energy_QMC_β2), marker=:x)
    # scatter!(Δ/2π, energy_ED[2,:], marker=:x)
end

@testset "1D Chain (9 atoms), β=0.5" begin
    β = 0.5
    EQ_MCS = 100
    MCS = 10_000
    M = 5000

    rng = MersenneTwister(1234)

    energy_QMC_β3 = []

    χ_squared = 0

    for ii in 1:Δ_step
        @show ii
        h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
        H = rydberg_QMC(h_ii)
        ts = BinaryThermalState(H, M)
        d = Diagnostics()
    
        [mc_step_beta!(rng, ts, H, β, d, eq=true) for i in 1:EQ_MCS] #equilibration phase
        
        ns = zeros(MCS)
    
        for i in 1:MCS # Monte Carlo Steps
            ns[i] = mc_step_beta!(rng, ts, H, β, d, eq=false)
        end

        # Binning analysis 
        energy(x) = -x / β + H.energy_shift
        BE = LogBinner(energy.(ns))
        τ_energy = tau(BE)
        @show τ_energy
        ratio = 2 * τ_energy + 1
        @show energy_binned = measurement(mean(BE), std_error(BE)*sqrt(ratio)) 
        append!(energy_QMC_β3, energy_binned)
        println()
        #append!(energy_QMC_β3, mean_and_stderr(x -> -x/β, ns) + H.energy_shift)

        @show χ_squared += abs2(value(energy_QMC_β3[ii]) - energy_ED[3, ii]) / abs2(uncertainty(energy_QMC_β3[ii]))
        # @test abs(stdscore(energy_QMC_β3[ii], energy_ED[2,ii])) < THRESHOLD_t
    end
    @test χ_squared < THRESHOLD_χ

    # scatter(Δ/2π, value.(energy_QMC_β3); yerror=uncertainty.(energy_QMC_β3), marker=:x)
    # scatter!(Δ/2π, energy_ED[3,:], marker=:x)
end

#########################

# Generate ED values for test with site-dependent parameters

β = 0.005
nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.48)

Ω_step = 9
Δ_step = 30
Ω = collect(LinRange(-2π * 0, 2π * 4, Ω_step))
Δ = LinRange(-2π * 9, 2π * 9, Δ_step)

energy_ED_sites = zeros(Δ_step)

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
    h_m = Matrix(mat(h_ii)) 
    energies, vecs = LinearAlgebra.eigen(h_m) 
    N_e = length(energies)

    w = exp.(-β .* (energies .- energies[1]))
    energy_ED_sites[ii] = sum(w .* energies) / sum(w)
end


@testset "1D Chain with site-dependent parameters" begin
    EQ_MCS = 30
    MCS = 100_000
    M = 50

    rng = MersenneTwister(3214)

    energy_QMC_sites = []

    χ_squared = 0

    for ii in 1:Δ_step
        @show ii
        h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
        H = rydberg_QMC(h_ii)
        ts = BinaryThermalState(H, M)
        d = Diagnostics()
    
        [mc_step_beta!(rng, ts, H, β, d, eq=true) for i in 1:EQ_MCS] #equilibration phase
    
        ns = zeros(MCS)
    
        for i in 1:MCS # Monte Carlo Steps
            ns[i] = mc_step_beta!(rng, ts, H, β, d, eq=false)
        end       

        # Binning analysis 
        energy(x) = -x / β + H.energy_shift
        BE = LogBinner(energy.(ns))
        τ_energy = tau(BE)
        ratio = 2 * τ_energy + 1
        energy_binned = measurement(mean(BE), std_error(BE)*sqrt(ratio)) 
        append!(energy_QMC_sites, energy_binned)

        χ_squared += abs2(value(energy_QMC_sites[ii]) - energy_ED_sites[ii]) / abs2(uncertainty(energy_QMC_sites[ii]))
        @show χ_squared
        # @test abs(stdscore(energy_QMC_β1[ii], energy_ED[1,ii])) < THRESHOLD_t       t-test for testing each QMC run individually
    end

    @test χ_squared < THRESHOLD_χ

    # scatter(Δ/2π, value.(energy_QMC_sites); yerror=uncertainty.(energy_QMC_sites), marker=:x)
    # scatter!(Δ/2π, energy_ED_sites, marker=:x)
end