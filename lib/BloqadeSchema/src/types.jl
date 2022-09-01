using GarishPrint


const ConstantParam = Union{Real,Nothing,Vector{<:Real}}
const DynamicParam = Union{Waveform{F,T} where {F,T<:Real},Vector{Waveform{F,T}} where {F,T<:Real} }


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

# Base.length(x::BloqadeSchema.RydbergDetuningLocal) = length(x.lattice_site_coefficients)

@option struct RydbergDetuning <: QuEraSchema
    global_value::RydbergDetuningGlobal
    local_value::Maybe{RydbergDetuningLocal}
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
    """
    Rabi Amplitude minimum value : 0 rad/s
    Rabi Amplitude maximum value : 25.0e6 rad/s
    Rabi Amplitude resolution : 400 rad/s
    Rabi Amplitude maximum slope : 2.5e14 (rad/s) / s
    
    Rabi Phase minimum value : -99 rad
    Rabi Phase maximum value : 99 rad
    Rabi Phase resolution : 0.5e-6 rad
    Rabi Phase maximum slope : 62.0e6 rad/s
    
    Detuning Amplitude minimum value : -125.0e6 rad/s
    Detuning Amplitude maximum value : 125.0e6 rad/s
    Detuning Amplitude resolution : 0.2 rad/s
    Detuning Amplitude maximum slope : 2.5e15 (rad/s) / s

    Local Detuning Scale Factor minimum value : 0
    Local Detuning Scale Factor maximum value : 1
    Local Detuning Scale Factor resolution : 0.01
    
    Time resolution : 1e-9 s
    Time step minimum : 1e-8 s
    """
    rabi_frequency_amplitude_maximum::Float64 = 0.0
    rabi_frequency_amplitude_minimum::Float64 = 25.0e6
    rabi_frequency_amplitude_resolution::Float64 = 400.0
    rabi_frequency_amplitude_max_slope::Float64 = 2.5e14

    rabi_frequency_phase_maximum::Float64 = -99.0
    rabi_frequency_phase_minimum::Float64 = 99.0
    rabi_frequency_phase_resolution::Float64 = 0.5e-6
    rabi_frequency_phase_max_slope::Float64 = 62.0e6

    rabi_detuning_maximum::Float64 = 125.0e6
    rabi_detuning_minimum::Float64 = -125.0e6
    rabi_detuning_resolution::Float64 = 0.2
    rabi_detuning_max_slope::Float64 = 2.5e15
    rabi_detuning_local_minimum::Float64 = 0.0
    rabi_detuning_local_maximum::Float64 = 1.0
    rabi_detuning_local_resolution::Float64 = 0.01

    rabi_time_resolution::Float64 = 1.0e-9
    rabi_time_min_step::Float64 = 1.0e-8
    rabi_time_maximum_value::Float64 = 4.0e-6
    
    waveform_tolerance::Float64 = 1.0e-3
    n_shots::Int = 1

    warn::Bool = false
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
