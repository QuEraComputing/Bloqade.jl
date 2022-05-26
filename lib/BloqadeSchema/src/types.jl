using GarishPrint

abstract type QuEraSchema end

@option struct Lattice <: QuEraSchema
    sites::Vector{Tuple{Float64,Float64}}
    filling::Vector{Int32}
end

@option struct RydbergRabiFrequencyAmplitudeGlobal <: QuEraSchema
    times::Vector{Float64}
    values::Vector{Float64}
end

@option struct RydbergRabiFrequencyAmplitude <: QuEraSchema
    global_value::RydbergRabiFrequencyAmplitudeGlobal
end

@option struct RydbergRabiFrequencyPhaseGlobal <: QuEraSchema
    times::Vector{Float64}
    values::Vector{Float64}
end

@option struct RydbergRabiFrequencyPhase <: QuEraSchema
    global_value::RydbergRabiFrequencyPhaseGlobal
end

@option struct RydbergDetuningGlobal <: QuEraSchema
    times::Vector{Float64}
    values::Vector{Float64}
end

@option struct RydbergDetuningLocal <: QuEraSchema
    times::Vector{Float64}
    values::Vector{Float64}
    lattice_site_coefficients::Vector{Float64}
end

@option struct RydbergDetuning <: QuEraSchema
    global_value::RydbergDetuningGlobal
    local_value::Maybe{Vector{RydbergDetuningLocal}}
end

@option struct RydbergHamiltonian <: QuEraSchema
    rabi_frequency_amplitude::RydbergRabiFrequencyAmplitude
    rabi_frequency_phase::RydbergRabiFrequencyPhase
    detuning::RydbergDetuning
end

@option struct EffectiveHamiltonian <: QuEraSchema
    rydberg::RydbergHamiltonian
end

@option struct TaskSpecification <: QuEraSchema
    nshots::Int
    lattice::Lattice
    effective_hamiltonian::EffectiveHamiltonian
end

@option struct SchemaConversionParams <: QuEraSchema
    rabi_frequency_amplitude_max_slope::Float64 = 0.1
    rabi_frequency_phase_max_slope::Float64 = 0.1
    rabi_detuning_max_slope::Float64 = 0.1
    n_shots::Int
end

@option struct ShotOutput <: QuEraSchema
    shot_status_code::Int32
    pre_sequence::Vector{Int32}
    post_sequence::Vector{Int32}
end

@option struct TaskOutput <: QuEraSchema
    task_status_code::Int
    shot_outputs::Vector{ShotOutput}
end

function Base.show(io::IO, ::MIME"text/plain", t::TaskSpecification)
    GarishPrint.pprint_struct(t)
end
