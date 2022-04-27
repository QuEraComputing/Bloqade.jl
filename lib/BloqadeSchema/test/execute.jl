using BloqadeSchema
using BloqadeExpr
using Test
using Configurations
using Yao

@testset "to_lattice with 1d chains" begin
    @test BloqadeSchema.to_lattice([1, 2, 3, 4, 5]) == BloqadeSchema.Lattice(;
        sites=[(1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
        filling=[1, 1, 1, 1, 1]
    )

    @test BloqadeSchema.to_lattice([]) == BloqadeSchema.Lattice(;
        sites=[],
        filling=[]
    )
end

@testset "to_lattice with 2d arrays" begin
    @test BloqadeSchema.to_lattice([(1, 3)]) == BloqadeSchema.Lattice(;
        sites=[(1, 3)],
        filling=[1]
    )

    @test BloqadeSchema.to_lattice([(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)]) == BloqadeSchema.Lattice(;
        sites=[(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)],
        filling=[1, 1, 1, 1, 1, 1]
    )
end

@testset "to_hamiltonian" begin
    Ω = BloqadeWaveforms.piecewise_constant(; clocks=[0, 2, 4, 6], values=[5, 3, 4, 6])
    Δ = BloqadeWaveforms.piecewise_linear(; clocks=[0.0, 0.6, 2.1, 2.2], values=[-10.1, -10.1, 10.1, 10.1])
    ϕ = BloqadeWaveforms.piecewise_linear(; clocks=[0, 5], values=[33, 0])

    result = BloqadeSchema.to_hamiltonian(; ϕ=ϕ, Ω=Ω, Δ=Δ)

    @test result == BloqadeSchema.EffectiveHamiltonian(; rydberg=BloqadeSchema.RydbergHamiltonian(;
        rabi_frequency_amplitude=BloqadeSchema.RydbergRabiFrequencyAmplitude(;
            global_value=BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal(;
                times=[0, 2, 4, 6],
                values=[5, 3, 4, 6]
            )
        ),
        rabi_frequency_phase=BloqadeSchema.RydbergRabiFrequencyPhase(;
            global_value=BloqadeSchema.RydbergRabiFrequencyPhaseGlobal(;
                times=[0, 5],
                values=[33, 0]
            )
        ),
        detuning=BloqadeSchema.RydbergDetuning(;
            global_value=BloqadeSchema.RydbergDetuningGlobal(;
                times=[0.0, 0.6, 2.1, 2.2],
                values=[-10.1, -10.1, 10.1, 10.1]
            )
        )
    ))
end

@testset "to_schema" begin
    # Ω = BloqadeExpr.SumOfX(; Ω=BloqadeWaveforms.piecewise_constant(; clocks=[0, 2, 4, 6], values=[5, 3, 4, 6]), nsites=5)
    Ω = BloqadeWaveforms.piecewise_constant(; clocks=[0, 2, 4, 6], values=[5, 3, 4, 6])
    # Δ = BloqadeExpr.SumOfN(; Δ=BloqadeWaveforms.piecewise_linear(; clocks=[0.0, 0.6, 2.1, 2.2], values=[-10.1, -10.1, 10.1, 10.1]), nsites=5)
    Δ = BloqadeWaveforms.piecewise_linear(; clocks=[0.0, 0.6, 2.1, 2.2], values=[-10.1, -10.1, 10.1, 10.1])
    # Φ = BloqadeExpr.SumOfXPhase(; ϕ=BloqadeWaveforms.piecewise_linear(; clocks=[0, 5], values=[33, 0]), nsites=5)
    ϕ = BloqadeWaveforms.piecewise_linear(; clocks=[0, 5], values=[33, 0])
    atoms = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)]
    # block = Yao.AbstractBlock(Ω, Δ, Φ, atoms)
    block = BloqadeExpr.rydberg_h(atoms; Δ=Δ, Ω=Ω, ϕ=ϕ)
    @test BloqadeSchema.to_schema(block) == TaskSpecification(;
        nshots=1,
        lattice=BloqadeSchema.Lattice(;
            sites=[(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)],
            filling=[1, 1, 1, 1, 1, 1]
        ),
        effective_hamiltonian=BloqadeSchema.EffectiveHamiltonian(; rydberg=BloqadeSchema.RydbergHamiltonian(;
            rabi_frequency_amplitude=BloqadeSchema.RydbergRabiFrequencyAmplitude(;
                global_value=BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal(;
                    times=[0, 2, 4, 6],
                    values=[5, 3, 4, 6]
                )
            ),
            rabi_frequency_phase=BloqadeSchema.RydbergRabiFrequencyPhase(;
                global_value=BloqadeSchema.RydbergRabiFrequencyPhaseGlobal(;
                    times=[0, 5],
                    values=[33, 0]
                )
            ),
            detuning=BloqadeSchema.RydbergDetuning(;
                global_value=BloqadeSchema.RydbergDetuningGlobal(;
                    times=[0.0, 0.6, 2.1, 2.2],
                    values=[-10.1, -10.1, 10.1, 10.1]
                )
            )
        ))
    )
end