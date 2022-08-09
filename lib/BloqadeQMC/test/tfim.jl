using Test, QMC
using Statistics
using Random
using RandomNumbers
using Measurements


expected_values = Dict{Tuple{Bool, Int, Float64}, Dict{String, Float64}}(
    (false, 10, 1.0) => Dict("M" => 6.451643863139102e-16,
                             "|M|" => 0.48188450771036156,
                             "M^2" => 0.32936397228034336,
                             "H" => -1.0774361711997833,
                             "C" => 3.085429815831006),

    (true, 10, 1.0) => Dict("M" => 1.4357874060005093e-16,
                            "|M|" => 0.5526982106525444,
                            "M^2" => 0.41362453160389673,
                            "H" => -1.1247786719410742,
                            "C" => 3.236282462798428),

    (true, 10, Inf64) => Dict("M" => -1.0480505352461478e-14,
                              "|M|" => 0.7072209347012945,
                              "M^2" => 0.5853102584019818,
                              "H" => -1.2784906442999324),

    (false, 10, Inf64) => Dict("M" => -1.5099033134902129e-15,
                               "|M|" => 0.5609773650317478,
                               "M^2" => 0.4115332272881885,
                               "H" => -1.2381489999654751)
)


# @testset "1D TFIM Ground State $N-sites, PBC=$PBC" for PBC in [true, false], N in [10]
#     rng = Xorshifts.Xoroshiro128Plus(1234)

#     bonds, Ns, Nb = lattice_bond_spins((10,), PBC)
#     H = TFIM(bonds, 1, Ns, Nb, 1.0, 1.0)
#     gs = BinaryGroundState(H, 1000)

#     MCS = 1_000_000
#     EQ_MCS = 100_000

#     mags = zeros(MCS)
#     ns = zeros(MCS)
#     nb = zeros(MCS)

#     [mc_step!(rng, gs, H) for i in 1:EQ_MCS]

#     for i in 1:MCS # Monte Carlo Production Steps
#         mc_step!(rng, gs, H) do lsize, gs, H
#             spin_prop = sample(H, gs)
#             ns[i] = num_single_site_diag(H, gs.operator_list)
#             nb[i] = num_two_site_diag(H, gs.operator_list)
#             mags[i] = magnetization(spin_prop)
#         end
#     end

#     abs_mag = mean_and_stderr(abs, mags)
#     mag_sqr = mean_and_stderr(abs2, mags)

#     energy = jackknife(ns) do n
#         if H.h != 0
#             (-H.h / n) + H.energy_shift / nspins(H)
#         else
#             H.energy_shift / nspins(H)
#         end
#     end

#     expected_vals = expected_values[(PBC, N, Inf64)]
#     @test stdscore(abs_mag, expected_vals["|M|"]) < THRESHOLD
#     @test stdscore(mag_sqr, expected_vals["M^2"]) < THRESHOLD
#     @test stdscore(energy, expected_vals["H"]) < THRESHOLD
# end


@testset "1D TFIM Thermal State $N-sites, PBC=$PBC, β=1.0" for PBC in [true, false], N in [10]
    rng = Xorshifts.Xoroshiro128Plus(1234)
    # rng = MersenneTwister(4321)

    bonds, Ns, Nb = lattice_bond_spins((10,), PBC)
    H = TFIM(bonds, 1, Ns, Nb, 1.0, 1.0)
    th = BinaryThermalState(H, 1000)
    beta = 1.0

    MCS = 1_000_000
    EQ_MCS = 100_000

    mags = zeros(MCS)
    ns = zeros(MCS)

    [mc_step_beta!(rng, th, H, beta) for i in 1:EQ_MCS]

    for i in 1:MCS # Monte Carlo Steps
        ns[i] = mc_step_beta!(rng, th, H, beta) do lsize, th, H
            mags[i] = magnetization(sample(H, th))
        end
    end
    abs_mag = mean_and_stderr(abs, mags)
    mag_sqr = mean_and_stderr(abs2, mags)

    energy = mean_and_stderr(x -> -x/beta, ns) + H.energy_shift
    energy /= nspins(H)

    heat_capacity = jackknife(ns .^ 2, ns) do nsqr, n
        nsqr - n^2 - n
    end

    expected_vals = expected_values[(PBC, N, 1.0)]
    @test abs(stdscore(abs_mag, expected_vals["|M|"])) < THRESHOLD
    @test abs(stdscore(mag_sqr, expected_vals["M^2"])) < THRESHOLD
    @test abs(stdscore(energy, expected_vals["H"])) < THRESHOLD
    @test abs(stdscore(heat_capacity, expected_vals["C"])) < THRESHOLD
end
