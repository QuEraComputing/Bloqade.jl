using Test, BloqadeQMC
using BloqadeQMC: Square, Kagome
using Statistics
using Random
using RandomNumbers
using Measurements
using BinningAnalysis

# ED values for groundstate simulation with (R_b, δ) = (1.7, 1.3)

expected_values = Dict{Tuple{String, Int, Int, Tuple{Bool, Bool}}, Dict{String, Float64}}(
    ("Square", 2, 2, (false, false)) => Dict("energy_density" => -1.0735253983492834, 
                                            "mag density" => -0.266628852869978),
    ("Kagome", 2, 2, (false, true)) => Dict("energy_density" => -1.0594971252029348, 
                                            "mag density" => -0.35758556716172124))

@testset "2D Square Rydberg Ground State, N = 2x2, PBC=(false, false)" begin
    rng = MersenneTwister(1234)

    R_b = 1.7
    delta = 3.3
    Ω = 1.0
    n1 = 2
    n2 = 2
    t = 1.0


    #lat = Kagome(t, n1, n2, (false, true); trunc=Inf)
    lat = Square(n1, n2, t, (false, false); trunc=Inf)
    H = Rydberg(lat, R_b, Ω, delta)
    gs = BinaryGroundState(H, 1000)         # Not exactly sure what the 1000 stands for
    d = Diagnostics()

    MCS = 1_000_000
    EQ_MCS = 100_000

    mags = zeros(MCS)
    ns = zeros(MCS)

    for i in 1:EQ_MCS
        mc_step!(rng, gs, H, d)
    end

    for i in 1:MCS # Monte Carlo Production Steps
        mc_step!(rng, gs, H, d) do lsize, gs, H
            spin_prop = sample(H, gs, 1)
            ns[i] = num_single_site_diag(H, gs.operator_list)   # ns[i] already returned by mc_step_beta!() in finite T case
            mags[i] = magnetization(spin_prop)   # Implementation same as for finite T
        end
    end

    abs_mag = mean_and_stderr(abs, mags)
    mag_sqr = mean_and_stderr(abs2, mags)

    println(abs_mag)
    println(mag_sqr)
    println(abs_mag / nspins(H))
    println(mag_sqr / nspins(H))

    # energy = jackknife(ns) do n
    #     if H.hx != 0
    #         (-H.hx * (1.0 / n)) + H.energy_shift / nspins(H)
    #     else
    #         H.energy_shift / nspins(H)
    #     end
    # end

end