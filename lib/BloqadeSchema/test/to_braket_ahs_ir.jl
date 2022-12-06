using BloqadeLattices
using BloqadeExpr
using BloqadeSchema
using BloqadeWaveforms
# using Braket
using DecFP

random_times  = collect(0:9)
random_values = random_coefficients = rand(-10:10e-5:100,10)
rand_sites = [site for site in zip(rand(0:1e-7:100,10), rand(0:1e-7:100,10))]
lattice = BloqadeSchema.Lattice(
    sites = rand_sites,
    filling = fill(1, length(rand_sites))
)

@testset "Lattice" begin

    setup = to_braket_ahs_ir(lattice)
    # filling should be equal to number of atom positions
    @test length(setup.ahs_register.filling) === length(setup.ahs_register.sites)
    # atom positions should be the same from the original 
    converted_sites = map(lattice.sites) do site
        Dec128.(string.([site...]))
    end
        
    @test converted_sites == setup.ahs_register.sites
end


local_detuning = BloqadeSchema.RydbergDetuningLocal(
    random_times,
    random_values,
    random_coefficients
)
@testset "RydbergDetuningLocal" begin
    
    # to_braket_ahs_ir(local_value::BloqadeSchema.RydbergDetuningLocal)

    physical_field = to_braket_ahs_ir(local_detuning)
    # check that values remain consistent
    
    @test Dec128.(string.(random_times)) == physical_field.time_series.times
    @test Dec128.(string.(random_values)) == physical_field.time_series.values
    @test Dec128.(string.(random_coefficients)) == physical_field.pattern

end

amplitude_global = BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal(
    random_times,
    random_values
)

global_phase = BloqadeSchema.RydbergRabiFrequencyPhaseGlobal(
    random_times,
    random_values
)
global_detuning = BloqadeSchema.RydbergDetuningGlobal(
    random_times,
    random_values
)
@testset "Globals" begin

    # to_braket_ahs_ir(global_values::Union{BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal,
    # BloqadeSchema.RydbergRabiFrequencyPhaseGlobal,
    # BloqadeSchema.RydbergDetuningGlobal})
    @testset "RydbergRabiFrequencyAmplitudeGlobal" begin

        physical_field = to_braket_ahs_ir(amplitude_global)
        @test Dec128.(string.(random_times)) == physical_field.time_series.times
        @test Dec128.(string.(random_values)) == physical_field.time_series.values
        @test "uniform" == physical_field.pattern
    end

    @testset "RydbergRabiFrequencyPhaseGlobal" begin
        physical_field = to_braket_ahs_ir(global_phase)
        @test Dec128.(string.(random_times)) == physical_field.time_series.times
        @test Dec128.(string.(random_values)) == physical_field.time_series.values
        @test "uniform" == physical_field.pattern
    end

    @testset "RydbergDetuningGlobal" begin
        physical_field = to_braket_ahs_ir(global_detuning)
        @test Dec128.(string.(random_times)) == physical_field.time_series.times
        @test Dec128.(string.(random_values)) == physical_field.time_series.values
        @test "uniform" == physical_field.pattern
    end

end

rydberg_amplitude = BloqadeSchema.RydbergRabiFrequencyAmplitude(amplitude_global)
rydberg_phase = BloqadeSchema.RydbergRabiFrequencyPhase(global_phase)
@testset "Unwrap Globals" begin

    # to_braket_ahs_ir(amplitude_or_phase::Union{BloqadeSchema.RydbergRabiFrequencyAmplitude,
    # BloqadeSchema.RydbergRabiFrequencyPhase})

    random_times  = collect(0:9)
    random_values = rand(-10:10e-5:100,10)

    @testset "RydbergRabiFrequencyAmplitude" begin

        physical_field = to_braket_ahs_ir(rydberg_amplitude)
        @test Dec128.(string.(random_times)) == physical_field.time_series.times
        @test Dec128.(string.(random_values)) == physical_field.time_series.values
        @test "uniform" == physical_field.pattern
    end

    @testset "RydbergRabiFrequencyPhase" begin

        physical_field = to_braket_ahs_ir(rydberg_phase)
        @test Dec128.(string.(random_times)) == physical_field.time_series.times
        @test Dec128.(string.(random_values)) == physical_field.time_series.values
        @test "uniform" == physical_field.pattern
    end
end

detunings_without_local = BloqadeSchema.RydbergDetuning(
                            global_detuning,
                            nothing
                        )

detunings_with_local = BloqadeSchema.RydbergDetuning(
                            global_detuning,
                            local_detuning
                        )
BloqadeSchema.RydbergDetuningLocal(
                            random_times, 
                            random_values,
                            random_coefficients
                        )
@testset "RydbergDetuning" begin
    # function to_braket_ahs_ir(detuning::BloqadeSchema.RydbergDetuning) 

    @testset "without local detuning" begin
        Δ, δ = to_braket_ahs_ir(detunings_without_local)
        # Δ should be 
        @test δ === nothing
        @test Dec128.(string.(random_times)) == Δ.time_series.times
        @test Dec128.(string.(random_values)) == Δ.time_series.values
        @test "uniform" == Δ.pattern

    end


    @testset "with local detuning" begin

        Δ, δ = to_braket_ahs_ir(detunings_with_local)

        @test Dec128.(string.(random_times)) == δ.time_series.times
        @test Dec128.(string.(random_values)) == δ.time_series.values
        @test Dec128.(string.(random_coefficients)) == δ.pattern

        @test Dec128.(string.(random_times)) == Δ.time_series.times
        @test Dec128.(string.(random_values)) == Δ.time_series.values
        @test "uniform" == Δ.pattern
    end
end

rydberg_hamiltonian = BloqadeSchema.RydbergHamiltonian(
                        rydberg_amplitude,
                        rydberg_phase,
                        detunings_with_local
                    )

effective_hamiltonian = BloqadeSchema.EffectiveHamiltonian(
                            rydberg_hamiltonian
                        )

@testset "EffectiveHamiltonian" begin 

    # function to_braket_ahs_ir(hamiltonian::BloqadeSchema.EffectiveHamiltonian)

    hamiltonian = to_braket_ahs_ir(effective_hamiltonian) 
    @test length(hamiltonian.drivingFields) === 1
    @test length(hamiltonian.shiftingFields) === 1

    h_amplitude = hamiltonian.drivingFields[1].amplitude
    h_phase     = hamiltonian.drivingFields[1].phase
    h_detuning  = hamiltonian.drivingFields[1].detuning

    @test Dec128.(string.(random_times))  == h_amplitude.time_series.times
    @test Dec128.(string.(random_values)) == h_amplitude.time_series.values
    @test "uniform" == h_amplitude.pattern

    @test Dec128.(string.(random_times))  == h_phase.time_series.times
    @test Dec128.(string.(random_values)) == h_phase.time_series.values
    @test "uniform" == h_phase.pattern

    @test Dec128.(string.(random_times))  == h_detuning.time_series.times
    @test Dec128.(string.(random_values)) == h_detuning.time_series.values
    @test "uniform" == h_detuning.pattern
    
end

@testset "RydbergHamiltonian" begin
    
    # function to_braket_ahs_ir(hamiltonian::BloqadeSchema.RydbergHamiltonian)

    hamiltonian = to_braket_ahs_ir(rydberg_hamiltonian) 
    @test length(hamiltonian.drivingFields) === 1
    @test length(hamiltonian.shiftingFields) === 1

    h_amplitude = hamiltonian.drivingFields[1].amplitude
    h_phase     = hamiltonian.drivingFields[1].phase
    h_detuning  = hamiltonian.drivingFields[1].detuning

    @test Dec128.(string.(random_times))  == h_amplitude.time_series.times
    @test Dec128.(string.(random_values)) == h_amplitude.time_series.values
    @test "uniform" == h_amplitude.pattern

    @test Dec128.(string.(random_times))  == h_phase.time_series.times
    @test Dec128.(string.(random_values)) == h_phase.time_series.values
    @test "uniform" == h_phase.pattern

    @test Dec128.(string.(random_times))  == h_detuning.time_series.times
    @test Dec128.(string.(random_values)) == h_detuning.time_series.values
    @test "uniform" == h_detuning.pattern
end

@testset "TaskSpecification" begin
    task_specification = BloqadeSchema.TaskSpecification(
                            1000,
                            lattice,
                            effective_hamiltonian
                        )

    ahs_program = to_braket_ahs_ir(task_specification)

    # check header
    @test "braket.ir.ahs.program" == ahs_program.braketSchemaHeader.name
    @test "1" == ahs_program.braketSchemaHeader.version

    # check positions
    setup = ahs_program.setup
    @test length(setup.ahs_register.filling) === length(setup.ahs_register.sites)
    converted_sites = map(lattice.sites) do site
        Dec128.(string.([site...]))
    end
    @test converted_sites == setup.ahs_register.sites

    # check hamiltonian
    hamiltonian = ahs_program.hamiltonian

    @test length(hamiltonian.drivingFields) === 1
    @test length(hamiltonian.shiftingFields) === 1

    h_amplitude = hamiltonian.drivingFields[1].amplitude
    h_phase     = hamiltonian.drivingFields[1].phase
    h_detuning  = hamiltonian.drivingFields[1].detuning

    @test Dec128.(string.(random_times))  == h_amplitude.time_series.times
    @test Dec128.(string.(random_values)) == h_amplitude.time_series.values
    @test "uniform" == h_amplitude.pattern

    @test Dec128.(string.(random_times))  == h_phase.time_series.times
    @test Dec128.(string.(random_values)) == h_phase.time_series.values
    @test "uniform" == h_phase.pattern

    @test Dec128.(string.(random_times))  == h_detuning.time_series.times
    @test Dec128.(string.(random_values)) == h_detuning.time_series.values
    @test "uniform" == h_detuning.pattern

end
