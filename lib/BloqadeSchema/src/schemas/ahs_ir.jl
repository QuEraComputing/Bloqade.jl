


abstract type QuEraSchema end

@option struct Lattice <: QuEraSchema
    sites::Vector{Tuple{Float64,Float64}}
    filling::Vector{Int32}
end

@option struct GlobalField <: QuEraSchema
    times::Vector{Float64}
    values::Vector{Float64}
end

@option struct LocalField <: QuEraSchema
    times::Vector{Float64}
    values::Vector{Float64}
    lattice_site_coefficients::Vector{Float64}
end

@option struct RabiFrequencyAmplitude <: QuEraSchema
    global_value::GlobalField
end

@option struct RabiFrequencyPhase <: QuEraSchema
    global_value::GlobalField
end

@option struct Detuning <: QuEraSchema
    global_value::GlobalField
    local_value::Maybe{LocalField}
end

@option struct RydbergHamiltonian <: QuEraSchema
    rabi_frequency_amplitude::RabiFrequencyAmplitude
    rabi_frequency_phase::RabiFrequencyPhase
    detuning::Detuning
end

@option struct EffectiveHamiltonian <: QuEraSchema
    rydberg::RydbergHamiltonian
end

"""
    struct QuEraTaskSpecification <: QuEraSchema

The schema representation of a task for the machine.

Is the output of [`to_schema`](@ref) and [`to_schema_no_validation`](@ref)
as well as input to [`execute(task::QuEraTaskSpecification)`](@ref).

# Fields
- `nshots::Int`: Number of shots (number of times hamiltonian is executed)
- `lattice::Lattice`: The Bravais lattice vectors and sites
- `effective_hamiltonian::EffectiveHamiltonian`: a `RydbergHamiltonian` instance 
wrapped inside an `EffectiveHamiltonian`
"""
@option struct QuEraTaskSpecification <: QuEraSchema
    nshots::Int
    lattice::Lattice
    effective_hamiltonian::EffectiveHamiltonian
end
