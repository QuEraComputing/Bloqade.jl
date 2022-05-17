using Configurations
using BloqadeSchema
using BloqadeExpr
using BloqadeWaveforms
using Test
using Yao
using OrderedCollections

@testset "to_lattice with 1d chains" begin
    @test BloqadeSchema.to_lattice([1, 2, 3, 4, 5]) ==
          BloqadeSchema.Lattice(; sites = [(1, 0), (2, 0), (3, 0), (4, 0), (5, 0)], filling = [1, 1, 1, 1, 1])

    @test BloqadeSchema.to_lattice([]) == BloqadeSchema.Lattice(; sites = [], filling = [])
end

@testset "to_lattice with 2d arrays" begin
    @test BloqadeSchema.to_lattice([(1, 3)]) == BloqadeSchema.Lattice(; sites = [(1, 3)], filling = [1])

    @test BloqadeSchema.to_lattice([(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)]) == BloqadeSchema.Lattice(;
        sites = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)],
        filling = [1, 1, 1, 1, 1, 1],
    )
end

@testset "to_hamiltonian" begin
    Ω = BloqadeWaveforms.piecewise_constant(; clocks = [0, 2, 4, 6, 6], values = [5, 3, 4, 6])
    Δ = BloqadeWaveforms.piecewise_linear(; clocks = [0.0, 0.6, 2.1, 2.2], values = [-10.1, -10.1, 10.1, 10.1])
    ϕ = BloqadeWaveforms.piecewise_linear(; clocks = [0, 5], values = [33, 0])

    result = BloqadeSchema.to_hamiltonian(Ω, ϕ, Δ, 10, 10, 10)

    @test result == BloqadeSchema.EffectiveHamiltonian(;
        rydberg = BloqadeSchema.RydbergHamiltonian(;
            rabi_frequency_amplitude = BloqadeSchema.RydbergRabiFrequencyAmplitude(;
                global_value = BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal(;
                    times = [0.0, 1.8, 2.0, 3.9, 4.0, 5.8, 6.0],
                    values = [5.0, 5.0, 3.0, 3.0, 4.0, 4.0, 6.0],
                ),
            ),
            rabi_frequency_phase = BloqadeSchema.RydbergRabiFrequencyPhase(;
                global_value = BloqadeSchema.RydbergRabiFrequencyPhaseGlobal(; times = [0, 5], values = [33, 0]),
            ),
            detuning = BloqadeSchema.RydbergDetuning(;
                global_value = BloqadeSchema.RydbergDetuningGlobal(;
                    times = [0.0, 0.6, 2.1, 2.2],
                    values = [-10.1, -10.1, 10.1, 10.1],
                ),
            ),
        ),
    )
end

@testset "to_schema" begin
    Ω = BloqadeWaveforms.piecewise_constant(; clocks = [0, 2, 4, 6, 7], values = [5, 3, 4, 6])
    Δ = BloqadeWaveforms.piecewise_linear(; clocks = [0.0, 0.6, 2.1, 2.2], values = [-10.1, -10.1, 10.1, 10.1])
    ϕ = BloqadeWaveforms.piecewise_linear(; clocks = [0, 5], values = [33, 0])
    atoms = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)]
    block = BloqadeExpr.rydberg_h(atoms; Δ = Δ, Ω = Ω, ϕ = ϕ)
    @test BloqadeSchema.to_schema(
        block;
        rabi_frequency_amplitude_max_slope = 10,
        rabi_frequency_phase_max_slope = 10,
        rabi_detuning_max_slope = 10,
        n_shots = 100,
    ) == TaskSpecification(;
        nshots = 100,
        lattice = BloqadeSchema.Lattice(;
            sites = [(0.0, 0.0), (1.0, 3.0), (4.0, 2.0), (6.0, 3.0), (0.0, 5.0), (2.0, 5.0)],
            filling = Int32[1, 1, 1, 1, 1, 1],
        ),
        effective_hamiltonian = BloqadeSchema.EffectiveHamiltonian(
            BloqadeSchema.RydbergHamiltonian(
                BloqadeSchema.RydbergRabiFrequencyAmplitude(
                    BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal(
                        [0.0, 1.8, 2.0, 3.9, 4.0, 5.8, 7.0],
                        [5.0, 5.0, 3.0, 3.0, 4.0, 4.0, 6.0],
                    ),
                ),
                BloqadeSchema.RydbergRabiFrequencyPhase(
                    BloqadeSchema.RydbergRabiFrequencyPhaseGlobal([0.0, 5.0], [33.0, 0.0]),
                ),
                BloqadeSchema.RydbergDetuning(
                    BloqadeSchema.RydbergDetuningGlobal([0.0, 0.6, 2.1, 2.2], [-10.1, -10.1, 10.1, 10.1]),
                    nothing,
                ),
            ),
        ),
    )
end

@testset "to_schema with all piecewise_linear waveforms" begin
    Ω = BloqadeWaveforms.piecewise_linear(; clocks = [0.0, 0.6, 2.1, 2.2], values = [-10.1, -10.1, 10.1, 10.1])
    Δ = BloqadeWaveforms.piecewise_linear(; clocks = [0.0, 0.6, 2.1, 2.2], values = [-10.1, -10.1, 10.1, 10.1])
    ϕ = BloqadeWaveforms.piecewise_linear(; clocks = [0, 5], values = [33, 0])
    atoms = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)]

    block = BloqadeExpr.rydberg_h(atoms; Δ = Δ, Ω = Ω, ϕ = ϕ)

    @test BloqadeSchema.to_schema(
        block;
        rabi_frequency_amplitude_max_slope = 10,
        rabi_frequency_phase_max_slope = 10,
        rabi_detuning_max_slope = 10,
        n_shots = 100,
    ) == TaskSpecification(;
        nshots = 100,
        lattice = BloqadeSchema.Lattice(;
            sites = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)],
            filling = [1, 1, 1, 1, 1, 1],
        ),
        effective_hamiltonian = BloqadeSchema.EffectiveHamiltonian(;
            rydberg = BloqadeSchema.RydbergHamiltonian(;
                rabi_frequency_amplitude = BloqadeSchema.RydbergRabiFrequencyAmplitude(;
                    global_value = BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal(;
                        times = [0.0, 0.6, 2.1, 2.2],
                        values = [-10.1, -10.1, 10.1, 10.1],
                    ),
                ),
                rabi_frequency_phase = BloqadeSchema.RydbergRabiFrequencyPhase(;
                    global_value = BloqadeSchema.RydbergRabiFrequencyPhaseGlobal(; times = [0, 5], values = [33, 0]),
                ),
                detuning = BloqadeSchema.RydbergDetuning(;
                    global_value = BloqadeSchema.RydbergDetuningGlobal(;
                        times = [0.0, 0.6, 2.1, 2.2],
                        values = [-10.1, -10.1, 10.1, 10.1],
                    ),
                ),
            ),
        ),
    )
end

@testset "to_schema without ϕ" begin
    Ω = BloqadeWaveforms.piecewise_linear(; clocks = [0.0, 0.6, 2.1, 2.2], values = [-10.1, -10.1, 10.1, 10.1])
    Δ = BloqadeWaveforms.piecewise_linear(; clocks = [0.0, 0.6, 2.1, 2.2], values = [-10.1, -10.1, 10.1, 10.1])
    atoms = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)]

    block = BloqadeExpr.rydberg_h(atoms; Δ = Δ, Ω = Ω)

    @test BloqadeSchema.to_schema(
        block;
        rabi_frequency_amplitude_max_slope = 10,
        rabi_frequency_phase_max_slope = 10,
        rabi_detuning_max_slope = 10,
        n_shots = 100,
    ) == TaskSpecification(;
        nshots = 100,
        lattice = BloqadeSchema.Lattice(;
            sites = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)],
            filling = [1, 1, 1, 1, 1, 1],
        ),
        effective_hamiltonian = BloqadeSchema.EffectiveHamiltonian(;
            rydberg = BloqadeSchema.RydbergHamiltonian(;
                rabi_frequency_amplitude = BloqadeSchema.RydbergRabiFrequencyAmplitude(;
                    global_value = BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal(;
                        times = [0.0, 0.6, 2.1, 2.2],
                        values = [-10.1, -10.1, 10.1, 10.1],
                    ),
                ),
                rabi_frequency_phase = BloqadeSchema.RydbergRabiFrequencyPhase(;
                    global_value = BloqadeSchema.RydbergRabiFrequencyPhaseGlobal(; times = [0], values = [0]),
                ),
                detuning = BloqadeSchema.RydbergDetuning(;
                    global_value = BloqadeSchema.RydbergDetuningGlobal(;
                        times = [0.0, 0.6, 2.1, 2.2],
                        values = [-10.1, -10.1, 10.1, 10.1],
                    ),
                ),
            ),
        ),
    )
end

@testset "get_piecewise_linear_times_and_clocks" begin
    Ω = BloqadeWaveforms.piecewise_constant(; clocks = [0, 2, 4, 6, 6], values = [5, 3, 4, 6])
    ϕ = BloqadeWaveforms.piecewise_constant(; clocks = [0, 2.1, 4.9, 6.3, 6.3], values = [5.7, 3, 4.9, 9])
    Ψ = BloqadeWaveforms.piecewise_constant(; clocks = [0, 2, 4, 6, 6], values = [5, 4, 4, 6])
    Δ = BloqadeWaveforms.piecewise_linear(; clocks = [0.0, 0.6, 2.1, 2.2], values = [-10.1, -10.1, 10.1, 10.1])

    @test BloqadeSchema.get_piecewise_linear_times_and_clocks(Δ, 1.0) ==
          ([0.0, 0.6, 2.1, 2.2], [-10.1, -10.1, 10.1, 10.1])
    @test BloqadeSchema.get_piecewise_linear_times_and_clocks(Ω, 10) ==
          ([0, 1.8, 2, 3.9, 4, 5.8, 6], [5, 5, 3, 3, 4, 4, 6])
    @test BloqadeSchema.get_piecewise_linear_times_and_clocks(Ω, 5) ==
          ([0, 1.6, 2, 3.8, 4, 5.6, 6], [5, 5, 3, 3, 4, 4, 6])
    @test BloqadeSchema.get_piecewise_linear_times_and_clocks(ϕ, 10) ==
          ([0.0, 1.83, 2.1, 4.71, 4.9, 5.89, 6.3], [5.7, 5.7, 3.0, 3.0, 4.9, 4.9, 9.0])
    @test BloqadeSchema.get_piecewise_linear_times_and_clocks(Ψ, 10) == ([0, 1.9, 2, 4, 5.8, 6], [5, 5, 4, 4, 4, 6])

    @test BloqadeSchema.get_piecewise_linear_times_and_clocks(nothing, 10) == ([0], [0])
end

@testset "to_json no local detuning" begin
    Ω = BloqadeWaveforms.piecewise_constant(; clocks = [0, 2, 4, 6, 6], values = [5, 3, 4, 6])
    Δ = BloqadeWaveforms.piecewise_linear(; clocks = [0.0, 0.6, 2.1, 2.2], values = [-10.1, -10.1, 10.1, 10.1])
    ϕ = BloqadeWaveforms.piecewise_linear(; clocks = [0, 5], values = [33, 0])
    atoms = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)]

    block = BloqadeExpr.rydberg_h(atoms; Δ = Δ, Ω = Ω, ϕ = ϕ)
    params = BloqadeSchema.SchemaConversionParams(
        rabi_frequency_amplitude_max_slope = 10,
        rabi_frequency_phase_max_slope = 10,
        rabi_detuning_max_slope = 10,
        n_shots = 100,
    )

    @test BloqadeSchema.to_json(block, params) ==
          "{\"nshots\":100,\"lattice\":{\"sites\":[[0.0,0.0],[1.0,3.0],[4.0,2.0],[6.0,3.0],[0.0,5.0],[2.0,5.0]],\"filling\":[1,1,1,1,1,1]},\"effective_hamiltonian\":{\"rydberg\":{\"rabi_frequency_amplitude\":{\"global\":{\"times\":[0.0,1.8,2.0,3.9,4.0,5.8,6.0],\"values\":[5.0,5.0,3.0,3.0,4.0,4.0,6.0]}},\"rabi_frequency_phase\":{\"global\":{\"times\":[0.0,5.0],\"values\":[33.0,0.0]}},\"detuning\":{\"global\":{\"times\":[0.0,0.6,2.1,2.2],\"values\":[-10.1,-10.1,10.1,10.1]}}}}}"
end

@testset "from_dict after to_dict" begin
    Ω = BloqadeWaveforms.piecewise_constant(; clocks = [0, 2, 4, 6, 7], values = [5, 3, 4, 6])
    Δ = BloqadeWaveforms.piecewise_linear(; clocks = [0.0, 0.6, 2.1, 2.2], values = [-10.1, -10.1, 10.1, 10.1])
    ϕ = BloqadeWaveforms.piecewise_linear(; clocks = [0, 5], values = [33, 0])
    atoms = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)]

    block = BloqadeExpr.rydberg_h(atoms; Δ = Δ, Ω = Ω, ϕ = ϕ)
    params = BloqadeSchema.SchemaConversionParams(
        rabi_frequency_amplitude_max_slope = 10,
        rabi_frequency_phase_max_slope = 10,
        rabi_detuning_max_slope = 10,
        n_shots = 100,
    )

    d = BloqadeSchema.to_dict(block, params)
    @test Configurations.from_dict(BloqadeSchema.TaskSpecification, d) == BloqadeSchema.to_schema(
        block;
        rabi_frequency_amplitude_max_slope = 10,
        rabi_frequency_phase_max_slope = 10,
        rabi_detuning_max_slope = 10,
        n_shots = 100,
    )
end
