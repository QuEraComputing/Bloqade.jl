using Test, BloqadeQMC
using Statistics
using Random
using RandomNumbers
using Measurements 
using BinningAnalysis
using BloqadeQMC: jackknife

expected_values = Dict{Tuple{Bool, Int, Float64}, Dict{String, Float64}}(
    (false, 10, Inf64) => Dict("M" => 0.9300787137837009,
                               "|M|" => 0.9300800988262218,
                               "M^2" => 0.879412342055611,
                               "H" => -2.0818751039065506),

    (false, 10, 1.0) => Dict("M" => 0.9152481757564336,
                             "|M|" => 0.9152717961127222,
                             "M^2" => 0.8564384998357667,
                             "H" => -2.051491413329706,
                             "C" => 1.5664917850898519),

    (true, 10, Inf64) => Dict("M" => 0.9424709532842519,
                              "|M|" => 0.9424716425300058,
                              "M^2" => 0.9001728915806667,
                              "H" => -2.1664826655166842),

    (true, 10, 1.0) => Dict("M" => 0.9339108974848992,
                            "|M|" => 0.9339228094722204,
                            "M^2" => 0.8868178763718234,
                            "H" => -2.1469676306665013,
                            "C" => 1.1680236801192336),

    (false, 5, Inf64) => Dict("M" => 0.9176830687768431,
                              "|M|" => 0.9183612642212718,
                              "M^2" => 0.8757671253668294,
                              "H" => -1.9972680527419864),

    (false, 5, 1.0) => Dict("M" => 0.8964725143596023,
                            "|M|" => 0.8991340140307242,
                            "M^2" => 0.8496226242573195,
                            "H" => -1.9559435941739651,
                            "C" => 0.984259674060965),

    (true, 5, Inf64) => Dict("M" => 0.9424692513081341,
                             "|M|" => 0.9427934100320888,
                             "M^2" => 0.9120933231955718,
                             "H" => -2.166482937646996),

    (true, 5, 1.0) => Dict("M" => 0.9337493597827439,
                           "|M|" => 0.9349728326000805,
                           "M^2" => 0.9014248358212859,
                           "H" => -2.1468390919494333,
                           "C" => 0.5884527794434007)
)


# @testset "1D LTFIM Ground State $N-sites, PBC=$PBC" for PBC in [true, false], N in [5, 10]
#     rng = Xorshifts.Xoroshiro128Plus(1234)

#     H = LTFIM((N,), 1.0, 1.0, 1.0, PBC)

#     gs = BinaryGroundState(H, 1000)

#     MCS = 1_000_000
#     EQ_MCS = 10_000

#     mags = zeros(MCS)
#     ns = zeros(MCS)

#     for i in 1:EQ_MCS
#         mc_step!(rng, gs, H)
#     end

#     for i in 1:MCS # Monte Carlo Production Steps
#         mc_step!(rng, gs, H) do lsize, gs, H
#             spin_prop = sample(H, gs)
#             ns[i] = num_single_site_diag(H, gs.operator_list)
#             mags[i] = magnetization(spin_prop)
#         end
#     end

#     abs_mag = mean_and_stderr(abs, mags)
#     mag_sqr = mean_and_stderr(abs2, mags)

#     energy = jackknife(ns) do n
#         if H.hx != 0
#             (-H.hx * (1.0 / n)) + H.energy_shift / nspins(H)
#         else
#             H.energy_shift / nspins(H)
#         end
#     end

#     expected_vals = expected_values[(PBC, N, Inf64)]
#     @test stdscore(abs_mag, expected_vals["|M|"]) < THRESHOLD
#     @test stdscore(mag_sqr, expected_vals["M^2"]) < THRESHOLD
#     @test stdscore(energy, expected_vals["H"]) < THRESHOLD
# end

THRESHOLD = 2.576  # 99% Two-sided CI of the t-distribution with infinite dofs

@testset "1D LTFIM Thermal State $N-sites, PBC=$PBC, β=1.0" for PBC in [true], N in [5]
    # rng = Xorshifts.Xoroshiro128Plus(1234)
    rng = MersenneTwister(2431)

    H = LTFIM((N,), 1.0, 1.0, 1.0, PBC)
    th = BinaryThermalState(H, 1000)
    beta = 1.0
    d = Diagnostics()

    MCS = 1_000_000
    EQ_MCS = 100_000

    mags = zeros(MCS)
    ns = zeros(MCS)

    [mc_step_beta!(rng, th, H, beta, d) for i in 1:EQ_MCS]

    for i in 1:MCS # Monte Carlo Steps
        ns[i] = mc_step_beta!(rng, th, H, beta, d) do lsize, th, H
            mags[i] = magnetization(sample(H, th, 1))
        end
    end

    # Binning Analysis for observables based on mags
    B = LogBinner(abs.(mags))
    println("Correlation time for abs_mag_binned samples")
    τ_abs = tau(B)
    @show τ_abs
    N_eff = MCS / (2 * τ_abs + 1)
    @show N_eff
    abs_mag_binned = measurement(mean(B), std_error(B)) 
    println()

    B2 = LogBinner(abs2.(mags))
    println("Correlation time for mag_sqr_binned samples")
    τ_abs2 = tau(B)
    @show τ_abs2
    N_eff2 = MCS / (2 * τ_abs2 + 1)
    @show N_eff2
    mag_sqr_binned = measurement(mean(B2), std_error(B2))
    println()

    # Binning analysis for observables based on ns
    energy(x) = -x / beta + H.energy_shift
    BE = LogBinner(energy.(ns) / nspins(H))
    println("Correlation time for energy_binned samples")
    τ_energy = tau(BE)
    @show τ_energy
    N_eff_E = MCS / (2 * τ_energy + 1)
    @show N_eff_E
    energy_density_binned = measurement(mean(BE), std_error(BE))
    println()

    # Unbinned calculations

    abs_mag = mean_and_stderr(abs, mags)
    @show abs_mag
    @show abs_mag_binned
    println()

    mag_sqr = mean_and_stderr(abs2, mags)
    @show mag_sqr
    @show mag_sqr_binned
    println()

    energy = mean_and_stderr(x -> -x/beta, ns) + H.energy_shift
    energy_density = energy / nspins(H)
    @show energy_density
    @show energy_density_binned
    println()


    heat_capacity = jackknife(ns .^ 2, ns) do nsqr, n
        nsqr - n^2 - n
    end

    @show heat_capacity

    expected_vals = expected_values[(PBC, N, 1.0)]
    # @test abs(stdscore(abs_mag, expected_vals["|M|"])) < THRESHOLD
    @test abs(stdscore(abs_mag_binned, expected_vals["|M|"])) < THRESHOLD
    # @test abs(stdscore(mag_sqr, expected_vals["M^2"])) < THRESHOLD
    @test abs(stdscore(mag_sqr_binned, expected_vals["M^2"])) < THRESHOLD
    # @test abs(stdscore(energy, expected_vals["H"])) < THRESHOLD
    @test abs(stdscore(energy_density_binned, expected_vals["H"])) < THRESHOLD
    @test abs(stdscore(heat_capacity, expected_vals["C"])) < THRESHOLD
end
