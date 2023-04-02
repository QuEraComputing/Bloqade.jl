using Test
using Configurations
using BloqadeSchema

task = QuEraTaskSpecification(;
    nshots = 239,
    lattice = BloqadeSchema.Lattice(; sites = [(0, 0), (5, 0), (0, 5), (5, 5)], filling = [1, 0, 1, 0]),
    effective_hamiltonian = BloqadeSchema.EffectiveHamiltonian(;
        rydberg = BloqadeSchema.RydbergHamiltonian(;
            rabi_frequency_amplitude = BloqadeSchema.RabiFrequencyAmplitude(;
                global_value = BloqadeSchema.GlobalField(;
                    times = [0, 5, 10],
                    values = [3, 4, 5],
                ),
            ),
            rabi_frequency_phase = BloqadeSchema.RabiFrequencyPhase(;
                global_value = BloqadeSchema.RabiFrequencyPhaseGlobal(; times = [0, 3, 6], values = [5, 4, 3]),
            ),
            detuning = BloqadeSchema.Detuning(;
                global_value = BloqadeSchema.GlobalField(; times = [0, 2, 4], values = [3, 4, 5]),
                local_value = [
                    BloqadeSchema.LocalField(;
                        times = [1, 3, 5],
                        values = [5, 4, 3],
                        lattice_site_coefficients = [1, 0, 1, 0],
                    ),
                    BloqadeSchema.LocalField(;
                        times = [5, 6, 7],
                        values = [3, 2, 1],
                        lattice_site_coefficients = [1, 1, 1, 1],
                    ),
                ],
            ),
        ),
    ),
)

to_dict(task)["effective_hamiltonian"]["rydberg"]["detuning"]["local"]
# to_toml("test.toml", task)
# from_toml(QuEraTaskSpecification, pkgdir(BloqadeSchema, "test", "test.toml"))
