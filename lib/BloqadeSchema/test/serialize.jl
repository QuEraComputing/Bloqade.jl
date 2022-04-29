using Test
using Configurations
using BloqadeSchema

task = TaskSpecification(;
    nshots=239,
    lattice=BloqadeSchema.Lattice(;
        sites=[(0, 0), (5, 0), (0, 5), (5, 5)],
        filling=[1, 0, 1, 0]
    ),
    effective_hamiltonian=BloqadeSchema.EffectiveHamiltonian(;
        rydberg=BloqadeSchema.RydbergHamiltonian(;
            rabi_frequency_amplitude=BloqadeSchema.RydbergRabiFrequencyAmplitude(;
                global_value=BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal(;
                    times=[0, 5, 10],
                    values=[3, 4, 5]
                )
            ),
            rabi_frequency_phase=BloqadeSchema.RydbergRabiFrequencyPhase(;
                global_value=BloqadeSchema.RydbergRabiFrequencyPhaseGlobal(;
                    times=[0, 3, 6],
                    values=[5, 4, 3]
                )
            ),
            detuning=BloqadeSchema.RydbergDetuning(;
                global_value=BloqadeSchema.RydbergDetuningGlobal(;
                    times=[0, 2, 4],
                    values=[3, 4, 5]
                ),
                local_value=[BloqadeSchema.RydbergDetuningLocal(;
                        times=[1, 3, 5],
                        values=[5, 4, 3],
                        lattice_site_coefficients=[1, 0, 1, 0]
                    ), BloqadeSchema.RydbergDetuningLocal(;
                        times=[5, 6, 7],
                        values=[3, 2, 1],
                        lattice_site_coefficients=[1, 1, 1, 1]
                    )]
            )
        )
    )
)

to_dict(task)["effective_hamiltonian"]["rydberg"]["detuning"]["local"]
# to_toml("test.toml", task)
# from_toml(TaskSpecification, pkgdir(BloqadeSchema, "test", "test.toml"))
