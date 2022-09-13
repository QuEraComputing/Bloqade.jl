using GarishPrint

const PiecewiseLinearWaveform = Waveform{BloqadeWaveforms.PiecewiseLinear{T,I},T} where {T<:Real,I}



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


# copying structure from device API 

# see https://github.com/QuEra-QCS/TaskManager/blob/main/api-impl/src/main/kotlin/com/queraqcs/api/impl/services/CapabilitiesStore.kt
"""
val capabilitiesStore = mapOf(
    "qpu1-mock" to VersionedCapabilities(
        version = "0.1",
        DeviceCapabilities(
            TaskCapabilities(
                numberQubitsMax = 256,
                numberShotsMin = 1,
                numberShotsMax = 10_000,
            ),
            LatticeCapabilities(
                LatticeAreaCapabilities(
                    width = 100e-6,
                    height = 100e-6,
                ),
                LatticeGeometryCapabilities(
                    spacingRadialMin = 4e-6,
                    spacingVerticalMin = 2.5e-6,
                    positionResolution = 0.1e-6,
                    numberSitesMax = 256,
                ),
            ),
            RydbergCapabilities(
                c6Coefficient = 5.420e-24,
                RydbergGlobalCapabilities(
                    rabiFrequencyMin = 0.0,
                    rabiFrequencyMax = 25.0e6,
                    rabiFrequencyResolution = 400.0,
                    rabiFrequencySlewRateMax = 2.5e14,
                    detuningMin = -125.0e6,
                    detuningMax = 125.0e6,
                    detuningResolution = 0.2,
                    detuningSlewRateMax = 2.5e15,
                    phaseMin = -99.0,
                    phaseMax = 99.0,
                    phaseResolution = 0.5e-6,
                    phaseSlewRateMax = 62.0e6,
                    timeMin = 0.0,
                    timeMax = 4e-6,
                    timeResolution = 1e-9,
                    timeDeltaMin = 10e-9,
                ),
                RydbergLocalCapabilities(
                    detuningMin = 0.0,
                    detuningMax = 125.0e6,
                    commonDetuningResolution = 2e3,
                    localDetuningResolution = 0.01,
                    detuningSlewRateMax = 1.25e15,
                    numberLocalDetuningSites = 256,
                    spacingRadialMin = 4e-6,
                    timeResolution = 1e-9,
                    timeDeltaMin = 10e-9,
                ),
            ),
        ),
    ),
)
"""

# make field mutable for user to change parameters

@option mutable struct LatticeAreaCapabilities <: QuEraSchema 
    width::Float64
    height::Float64
end

@option mutable struct LatticeGeometryCapabilities <: QuEraSchema 
    spacingRadialMin::Float64
    spacingVerticalMin::Float64
    positionResolution::Float64
    numberSitesMax::Int
end

@option struct LatticeCapabilities <: QuEraSchema 
    area::LatticeAreaCapabilities
    geometry::LatticeGeometryCapabilities
end

# non snake case is from API
@option mutable struct RydbergGlobalCapabilities <: QuEraSchema 
    rabiFrequencyMin::Float64
    rabiFrequencyMax::Float64
    rabiFrequencyResolution::Float64
    rabiFrequencySlewRateMax::Float64
    detuningMin::Float64
    detuningMax::Float64
    detuningResolution::Float64
    detuningSlewRateMax::Float64
    phaseMin::Float64
    phaseMax::Float64
    phaseResolution::Float64
    phaseSlewRateMax::Float64
    timeMin::Float64
    timeMax::Float64
    timeResolution::Float64
    timeDeltaMin::Float64
end

@option mutable struct RydbergLocalCapabilities <: QuEraSchema 
    detuningMin::Float64
    detuningMax::Float64
    commonDetuningResolution::Float64
    localDetuningResolution::Float64
    detuningSlewRateMax::Float64
    numberLocalDetuningSites::Int
    spacingRadialMin::Float64
    timeResolution::Float64
    timeDeltaMin::Float64
end

@option struct RydbergCapabilities <: QuEraSchema 
    c6Coefficient::Float64
    global_value::RydbergGlobalCapabilities
    local_value::RydbergLocalCapabilities
end

@option mutable struct TaskCapabilities <: QuEraSchema 
    numberQubitsMax::Int
    numberShotsMin::Int 
    numberShotsMax::Int
end

@option struct DeviceCapabilities <: QuEraSchema
    task::TaskCapabilities
    lattice::LatticeCapabilities
    rydberg::RydbergCapabilities
end


# manually convert to default units
get_device_capabilities() = DeviceCapabilities(
    task=TaskCapabilities(
        numberQubitsMax = 256,
        numberShotsMin = 1,
        numberShotsMax = 10000,
    ),
    lattice=LatticeCapabilities(
        area=LatticeAreaCapabilities(
            width = convert_units(100e-6,m,μm),
            height = convert_units(100e-6,m,μm)
        ),
        geometry=LatticeGeometryCapabilities(
            spacingRadialMin = convert_units(4e-6,m,μm),
            spacingVerticalMin = convert_units(2.5e-6,m,μm),
            positionResolution = convert_units(0.1e-6,m,μm),
            numberSitesMax = 256,
        )
    ),
    rydberg=RydbergCapabilities(
        c6Coefficient = convert_units(5.420e-24,rad*m^6/s,rad*μm^6/μs),
        global_value=RydbergGlobalCapabilities(
            rabiFrequencyMin = convert_units(0.0,rad/s,rad*MHz),
            rabiFrequencyMax = convert_units(25.0e6,rad/s,rad*MHz),
            rabiFrequencyResolution = convert_units(400.0,rad/s,rad*MHz),
            rabiFrequencySlewRateMax = convert_units(2.5e14,rad/s^2,rad*MHz/μs),
            detuningMin = convert_units(-125.0e6,rad/s,rad*MHz),
            detuningMax = convert_units(125.0e6,rad/s,rad*MHz),
            detuningResolution = convert_units(0.2,rad/s,rad*MHz),
            detuningSlewRateMax = convert_units(2.5e15,rad/s^2,rad*MHz/μs),
            phaseMin = convert_units(-99.0,rad,rad),
            phaseMax = convert_units(99.0,rad,rad),
            phaseResolution = convert_units(0.5e-6,rad,rad),
            phaseSlewRateMax = convert_units(62.0e6,rad/s,rad/μs),
            timeMin = convert_units(0.0,s,μs),
            timeMax = convert_units(4e-6,s,μs),
            timeResolution = convert_units(1e-9,s,μs),
            timeDeltaMin = convert_units(10e-9,s,μs)
        ),
        local_value=RydbergLocalCapabilities(
            detuningMin = convert_units(0.0,rad/s,rad*MHz),
            detuningMax = convert_units(125.0e6,rad/s,rad*MHz),
            commonDetuningResolution = convert_units(2e3,rad/s,rad*MHz),
            localDetuningResolution = 0.01,
            detuningSlewRateMax = convert_units(1.25e15,rad/s^2,rad*MHz/μs),
            numberLocalDetuningSites = 256,
            spacingRadialMin = convert_units(4e-6,m,μm),
            timeResolution = convert_units(1e-9,s,μs),
            timeDeltaMin = convert_units(10e-9,s,μs)
        )
    )
)
# leave as SI units, needed for rounding purposes
# get_device_capabilities_SI() = DeviceCapabilities(
#     task=TaskCapabilities(
#         numberQubitsMax = 256,
#         numberShotsMin = 1,
#         numberShotsMax = 10000,
#     ),
#     lattice=LatticeCapabilities(
#         area=LatticeAreaCapabilities(
#             width = 100e-6,
#             height = 100e-6
#         ),
#         geometry=LatticeGeometryCapabilities(
#             spacingRadialMin = 4e-6,
#             spacingVerticalMin = 2.5e-6,
#             positionResolution = 0.1e-6,
#             numberSitesMax = 256,
#         )
#     ),
#     rydberg=RydbergCapabilities(
#         c6Coefficient = 5.420e-24,
#         global_value=RydbergGlobalCapabilities(
#             rabiFrequencyMin = 0.0,
#             rabiFrequencyMax = 25.0e6,
#             rabiFrequencyResolution = 400.0,
#             rabiFrequencySlewRateMax = 2.5e14,
#             detuningMin = -125.0e6,
#             detuningMax = 125.0e6,
#             detuningResolution = 0.2,
#             detuningSlewRateMax = 2.5e15,
#             phaseMin = -99.0,
#             phaseMax = 99.0,
#             phaseResolution = 0.5e-6,
#             phaseSlewRateMax = 62.0e6,
#             timeMin = 0.0,
#             timeMax = 4e-6,
#             timeResolution = 1e-9,
#             timeDeltaMin = 10e-9
#         ),
#         local_value=RydbergLocalCapabilities(
#             detuningMin = 0.0,
#             detuningMax = 125.0e6,
#             commonDetuningResolution = 2e3,
#             localDetuningResolution = 0.01,
#             detuningSlewRateMax = 1.25e15,
#             numberLocalDetuningSites = 256,
#             spacingRadialMin = 4e-6,
#             timeResolution = 1e-9,
#             timeDeltaMin = 10e-9
#         )
#     )
# )


function get_rydberg_capabilities(;device_capabilities::DeviceCapabilities=get_device_capabilities())
    return (
        Ω=(
            min_time_step = device_capabilities.rydberg.global_value.timeDeltaMin,
            time_resolution = device_capabilities.rydberg.global_value.timeResolution,
            max_time = device_capabilities.rydberg.global_value.timeMax,
            max_value = device_capabilities.rydberg.global_value.rabiFrequencyMax,
            min_value = device_capabilities.rydberg.global_value.rabiFrequencyMin,
            max_slope = device_capabilities.rydberg.global_value.rabiFrequencySlewRateMax,
            value_resolution = device_capabilities.rydberg.global_value.rabiFrequencyResolution
        ),
        ϕ = (
            min_time_step = device_capabilities.rydberg.global_value.timeDeltaMin,
            time_resolution = device_capabilities.rydberg.global_value.timeResolution,
            max_time = device_capabilities.rydberg.global_value.timeMax,
            max_value = device_capabilities.rydberg.global_value.phaseMax,
            min_value = device_capabilities.rydberg.global_value.phaseMin,
            max_slope = device_capabilities.rydberg.global_value.phaseSlewRateMax,
            value_resolution = device_capabilities.rydberg.global_value.phaseResolution
        ),
        Δ = (
            min_time_step = device_capabilities.rydberg.global_value.timeDeltaMin,
            time_resolution = device_capabilities.rydberg.global_value.timeResolution,
            max_time = device_capabilities.rydberg.global_value.timeMax,
            max_value = device_capabilities.rydberg.global_value.detuningMax,
            min_value = device_capabilities.rydberg.global_value.detuningMin,
            max_slope = device_capabilities.rydberg.global_value.detuningSlewRateMax,
            value_resolution = device_capabilities.rydberg.global_value.detuningResolution
        ),
        δ = (
            min_time_step = device_capabilities.rydberg.local_value.timeDeltaMin,
            time_resolution = device_capabilities.rydberg.local_value.timeResolution,
            max_time = device_capabilities.rydberg.local_value.timeMax,
            max_value = device_capabilities.rydberg.local_value.detuningMax,
            min_value = device_capabilities.rydberg.local_value.detuningMin,
            max_slope = device_capabilities.rydberg.local_value.detuningSlewRateMax,
            value_resolution = device_capabilities.rydberg.local_value.commonDetuningResolution,
            local_mask_resolution = device_capabilities.rydberg.local_value.localDetuningResolution
        )
    
    )
end


@option struct SchemaTranslationParams <: QuEraSchema
    n_shots::Int = 1
    warn::Bool=true
    transform_info::Bool=true
    device_capabilities::DeviceCapabilities = get_device_capabilities()
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
