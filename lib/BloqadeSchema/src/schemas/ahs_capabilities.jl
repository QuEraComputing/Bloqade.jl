

abstract type QuEraCapabilities end

# copying structure from device API 

# make field mutable for user to change parameters

@option mutable struct LatticeAreaCapabilities <: QuEraCapabilities 
    width::Float64
    height::Float64
end

@option mutable struct LatticeGeometryCapabilities <: QuEraCapabilities 
    spacing_radial_min::Float64
    spacing_vertical_min::Float64
    position_resolution::Float64
    number_sites_max::Int
end

@option mutable struct LatticeCapabilities <: QuEraCapabilities 
    area::LatticeAreaCapabilities
    geometry::LatticeGeometryCapabilities
    number_qubits_max::Int
end

# non snake case is from API
@option mutable struct RydbergGlobalCapabilities <: QuEraCapabilities 
    rabi_frequency_min::Float64
    rabi_frequency_max::Float64
    rabi_frequency_resolution::Float64
    rabi_frequency_slew_rate_max::Float64
    detuning_min::Float64
    detuning_max::Float64
    detuning_resolution::Float64
    detuning_slew_rate_max::Float64
    phase_min::Float64
    phase_max::Float64
    phase_resolution::Float64
    time_min::Float64
    time_max::Float64
    time_resolution::Float64
    time_delta_min::Float64
end

@option mutable struct RydbergLocalCapabilities <: QuEraCapabilities 
    detuning_min::Float64
    detuning_max::Float64
    common_detuning_resolution::Float64
    local_detuning_resolution::Float64
    detuning_slew_rate_max::Float64
    number_local_detuning_sites::Int
    spacing_radial_min::Float64
    time_resolution::Float64
    time_delta_min::Float64
end

@option struct RydbergCapabilities <: QuEraCapabilities 
    c6_coefficient::Float64
    global_value::RydbergGlobalCapabilities
    local_value::Maybe{RydbergLocalCapabilities}
end

@option mutable struct TaskCapabilities <: QuEraCapabilities 
    number_shots_min::Int 
    number_shots_max::Int
end

@option struct DeviceCapabilities <: QuEraCapabilities
    task::TaskCapabilities
    lattice::LatticeCapabilities
    rydberg::RydbergCapabilities
end
