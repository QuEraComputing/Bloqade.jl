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
    rabi_frequency_amplitude_maximum::Quantity = Quantity(0,rad/s)
    rabi_frequency_amplitude_minimum::Quantity = Quantity(25.0e6,rad/s)
    rabi_frequency_amplitude_resolution::Quantity = Quantity(400.0,rad/s)
    rabi_frequency_amplitude_max_slope::Quantity = Quantity(2.5e14,rad/s^2)

    rabi_frequency_phase_maximum::Quantity = Quantity(-99.0,rad)
    rabi_frequency_phase_minimum::Quantity = Quantity(99.0,rad)
    rabi_frequency_phase_resolution::Quantity = Quantity(0.5e-6,rad)
    rabi_frequency_phase_max_slope::Quantity = Quantity(62.0e6,rad/s)

    rabi_detuning_maximum::Quantity = Quantity(125.0e6,rad/s)
    rabi_detuning_minimum::Quantity = Quantity(-125.0e6,rad/s)
    rabi_detuning_resolution::Quantity = Quantity(0.2,rad/s)
    rabi_detuning_max_slope::Quantity = Quantity(2.5e15,rad/s^2)
    rabi_detuning_local_minimum::Quantity = Quantity(0.0,NoUnits)
    rabi_detuning_local_maximum::Quantity = Quantity(1.0,NoUnits)
    rabi_detuning_local_resolution::Quantity = Quantity(0.01,NoUnits)

    rabi_time_resolution::Quantity = Quantity(1.0e-9,s)
    rabi_time_min_step::Quantity = Quantity(1.0e-8,s)
    rabi_time_maximum_values::Quantity = Quantity(4.0e-6,s)
    
    waveform_tolerance::Float64 = 1.0e-3
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
