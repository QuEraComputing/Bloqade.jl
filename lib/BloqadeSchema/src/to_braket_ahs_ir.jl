"""
    function to_braket_ahs_ir(local_value::BloqadeSchema.RydbergDetuningLocal)

Converts `local_value` to `Braket.IR.PhysicalField`
"""
function to_braket_ahs_ir(local_value::BloqadeSchema.RydbergDetuningLocal)

    # should also be a PhysicalField
    time_series = Braket.IR.TimeSeries(
        Dec128.(string.(local_value.values)),
        Dec128.(string.(local_value.times)),
    )

    Braket.IR.PhysicalField(
        time_series,
        Dec128.(string.(local_value.lattice_site_coefficients))
    )
end

"""
    function to_braket_ahs_ir(global_values::Union{BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal,
    BloqadeSchema.RydbergRabiFrequencyPhaseGlobal,
    BloqadeSchema.RydbergDetuningGlobal})

Converts `global_values` to `Braket.IR.PhysicalField`
"""
function to_braket_ahs_ir(global_values::Union{BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal,
                                               BloqadeSchema.RydbergRabiFrequencyPhaseGlobal,
                                               BloqadeSchema.RydbergDetuningGlobal})

    time_series = Braket.IR.TimeSeries(
        Dec128.(string.(global_values.values)),
        Dec128.(string.(global_values.times)),
    )
    # Create PhysicalField
    Braket.IR.PhysicalField(
        time_series,
        "uniform"
    )
end

"""
    function to_braket_ahs_ir(amplitude_or_phase::Union{BloqadeSchema.RydbergRabiFrequencyAmplitude,
    BloqadeSchema.RydbergRabiFrequencyPhase})
                                           
Unwraps `amplitude_or_phase` to extract their `global_value` fields which are 
immediately evaluated by `function to_braket_ahs_ir(global_values::Union{BloqadeSchema.RydbergRabiFrequencyAmplitudeGlobal,
BloqadeSchema.RydbergRabiFrequencyPhaseGlobal,
BloqadeSchema.RydbergDetuningGlobal})`
"""
to_braket_ahs_ir(amplitude_or_phase::Union{BloqadeSchema.RydbergRabiFrequencyAmplitude,
                                           BloqadeSchema.RydbergRabiFrequencyPhase}) = to_braket_ahs_ir(amplitude_or_phase.global_value)

"""
    function to_braket_ahs_ir(detuning::BloqadeSchema.RydbergDetuning)

Converts `detuning` to `Braket.IR.PhysicalField`'s for the `global_value`
and `local_value` fields of `detuning`. Returns the converted `global_value` first 
followed by the `local_value`. 
``
"""
function to_braket_ahs_ir(detuning::BloqadeSchema.RydbergDetuning)
    Δ = to_braket_ahs_ir(detuning.global_value)
    # local_value is Maybe{RydbergDetuningLocal} so chance that 
    # it's `nothing` or the actual value
    δ = nothing
    if detuning.local_value !== nothing
        δ = to_braket_ahs_ir(detuning.local_value)
    end

    return Δ, δ
end

"""
    function to_braket_ahs_ir(rydberg_hamiltonian::BloqadeSchema.RydbergHamiltonian)

Converts `rydberg_hamiltonian` to `Braket.IR.Hamiltonian`
"""
function to_braket_ahs_ir(rydberg_hamiltonian::BloqadeSchema.RydbergHamiltonian)
    # PhysicalFields
    Ω = to_braket_ahs_ir(rydberg_hamiltonian.rabi_frequency_amplitude)
    ϕ = to_braket_ahs_ir(rydberg_hamiltonian.rabi_frequency_phase)
    Δ, δ = to_braket_ahs_ir(rydberg_hamiltonian.detuning)

    driving_field = [Braket.IR.DrivingField(Ω,ϕ,Δ)]
    shifting_field = Braket.IR.ShiftingField[]
    
    # check if local detuning is present
    if δ !== nothing 
        shifting_field = [Braket.IR.ShiftingField(δ)]
    end

    Braket.IR.Hamiltonian(driving_field,shifting_field)

end

"""
    function to_braket_ahs_ir(lattice::BloqadeSchema.Lattice)

Converts `lattice` to a `Braket.IR.Setup` instance.
"""
function to_braket_ahs_ir(lattice::BloqadeSchema.Lattice)
    # Convert lattice.sites to proper time for AtomArrangement, 
    # go from Vector{Tuple{Float64, Float64}} to Vector{Vector{Dec128}}
    type_converted_sites = map(collect.(lattice.sites)) do position
        map(position) do x
            Dec128(string(x))
        end
    end
    # type_converted_sites = convert(Vector{Vector{Dec128}}, collect.(lattice.sites)) 
    # create AtomArrangement for Braket IR, goes into 
    # a Setup object
    atom_arrangement = Braket.IR.AtomArrangement(
                            type_converted_sites,
                            lattice.filling
                        )
    # AtomArrangement assigned to ahs_register field
    Braket.IR.Setup(atom_arrangement)
end

"""
    function to_braket_ahs_ir(hamiltonian::BloqadeSchema.EffectiveHamiltonian)

Unwraps the `BloqadeSchema.RydbergHamiltonian` contained inside `hamiltonian`
and immediately evaluates it using
 `to_braket_ahs_ir(rydberg_hamiltonian::BloqadeSchema.RydbergHamiltonian)`
"""
to_braket_ahs_ir(hamiltonian::BloqadeSchema.EffectiveHamiltonian) = to_braket_ahs_ir(hamiltonian.rydberg::BloqadeSchema.RydbergHamiltonian)

# Public API, *THE* main function users should use
"""
   function to_braket_ahs_ir(bloqade_task::BloqadeSchema.TaskSpecification)

Converts a `bloqade_task` into a `Braket.IR.AHSProgram`
capable of being submitted to AWS Braket for execution on a QPU.

NOTE: `BloqadeSchema.TaskSpecification` contains a field `nshots` which
does not have a corresponding field in `Braket.IR.AHSProgram`. This value must be
fed as a keyword argument to a `Braket.AwsDevice` instance separately.
"""
function to_braket_ahs_ir(bloqade_task::BloqadeSchema.TaskSpecification)
    braket_schema_header = Braket.IR.braketSchemaHeader("braket.ir.ahs.program",
                                                        "1")
    setup = to_braket_ahs_ir(bloqade_task.lattice)
    hamiltonian = to_braket_ahs_ir(bloqade_task.effective_hamiltonian)

    Braket.IR.AHSProgram(braket_schema_header, setup, hamiltonian)
end