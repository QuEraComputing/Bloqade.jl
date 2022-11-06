


function parse_unit(curr::Module, unit::AbstractString)
    ex = Meta.parse(unit)
    unitmods = [Unitful]
    for m in Unitful.unitmodules
        # Find registered unit extension modules which are also loaded by
        # __module__ (required so that precompilation will work).
        if isdefined(curr, nameof(m)) && getfield(curr, nameof(m)) === m
            push!(unitmods, m)
        end
    end
    return Unitful.lookup_units(unitmods, ex)
end

function convert_units_recursive(input_value,units)

    if typeof(input_value) <: AbstractDict
        # make a copy
        new_values = Dict(input_value)
        # for each k,v in input_dict, 
        # find the unit in the units dict and convert v to that unit
        for (key,value) in input_value
            new_values[key] = convert_units_recursive(value,units[key])
        end
        return new_values
    else
        units == "NoUnits" && return input_value
        
        # Inside the units dictionary, 
        # replace instances of μm with m, 
        # THEN replace instances of s with μs
        to_units = replace(
                replace(units,
                "m"=>"μm"),
            "s"=>"μs")

        from = Base.eval(@__MODULE__,parse_unit(@__MODULE__,units))
        to   = Base.eval(@__MODULE__,parse_unit(@__MODULE__,to_units))

        return convert_units(input_value,from,to)
    end
end

function get_device_capabilities(capabilities_file=nothing)
    units_file = Base.Filesystem.joinpath(pkgdir(@__MODULE__),"config","capabilities-qpu1-mock-units.json")
    capabilities_json_units = JSON.parse(JSON.open(units_file))["capabilities"]

    capabilities_json_SI = if isnothing(capabilities_file)
        capabilities_file = Base.Filesystem.joinpath(pkgdir(@__MODULE__),"config","capabilities-qpu1-mock.json")
        JSON.parse(JSON.open(capabilities_file))["capabilities"]
    else
        JSON.parse(JSON.open(capabilities_file))
    end

    capabilities_json = convert_units_recursive(capabilities_json_SI,capabilities_json_units)

    return Configurations.from_dict(DeviceCapabilities,capabilities_json)
end


# leave as SI units, needed for rounding purposes

function get_device_capabilities_SI(capabilities_file=nothing)

    capabilities_json = if isnothing(capabilities_file)
        capabilities_file = Base.Filesystem.joinpath(pkgdir(@__MODULE__),"config","capabilities-qpu1-mock.json")
        JSON.parse(JSON.open(capabilities_file))["capabilities"]
    else
        JSON.parse(JSON.open(capabilities_file))
    end

    return Configurations.from_dict(DeviceCapabilities,capabilities_json)
end



function get_rydberg_capabilities(;device_capabilities::DeviceCapabilities=get_device_capabilities())
    return if isnothing(device_capabilities.rydberg.local_value)
        (
            Ω=(
                min_time_step = device_capabilities.rydberg.global_value.time_delta_min,
                time_resolution = device_capabilities.rydberg.global_value.time_resolution,
                max_time = device_capabilities.rydberg.global_value.time_max,
                max_value = device_capabilities.rydberg.global_value.rabi_frequency_max,
                min_value = device_capabilities.rydberg.global_value.rabi_frequency_min,
                max_slope = device_capabilities.rydberg.global_value.rabi_frequency_slew_rate_max,
                value_resolution = device_capabilities.rydberg.global_value.rabi_frequency_resolution
            ),
            ϕ = (
                min_time_step = device_capabilities.rydberg.global_value.time_delta_min,
                time_resolution = device_capabilities.rydberg.global_value.time_resolution,
                max_time = device_capabilities.rydberg.global_value.time_max,
                max_value = device_capabilities.rydberg.global_value.phase_max,
                min_value = device_capabilities.rydberg.global_value.phase_min,
                value_resolution = device_capabilities.rydberg.global_value.phase_resolution
            ),
            Δ = (
                min_time_step = device_capabilities.rydberg.global_value.time_delta_min,
                time_resolution = device_capabilities.rydberg.global_value.time_resolution,
                max_time = device_capabilities.rydberg.global_value.time_max,
                max_value = device_capabilities.rydberg.global_value.detuning_max,
                min_value = device_capabilities.rydberg.global_value.detuning_min,
                max_slope = device_capabilities.rydberg.global_value.detuning_slew_rate_max,
                value_resolution = device_capabilities.rydberg.global_value.detuning_resolution
            ),
            δ = nothing
        )
    else
        (
            Ω=(
                min_time_step = device_capabilities.rydberg.global_value.time_delta_min,
                time_resolution = device_capabilities.rydberg.global_value.time_resolution,
                max_time = device_capabilities.rydberg.global_value.time_max,
                max_value = device_capabilities.rydberg.global_value.rabi_frequency_max,
                min_value = device_capabilities.rydberg.global_value.rabi_frequency_min,
                max_slope = device_capabilities.rydberg.global_value.rabi_frequency_slew_rate_max,
                value_resolution = device_capabilities.rydberg.global_value.rabi_frequency_resolution
            ),
            ϕ = (
                min_time_step = device_capabilities.rydberg.global_value.time_delta_min,
                time_resolution = device_capabilities.rydberg.global_value.time_resolution,
                max_time = device_capabilities.rydberg.global_value.time_max,
                max_value = device_capabilities.rydberg.global_value.phase_max,
                min_value = device_capabilities.rydberg.global_value.phase_min,
                value_resolution = device_capabilities.rydberg.global_value.phase_resolution
            ),
            Δ = (
                min_time_step = device_capabilities.rydberg.global_value.time_delta_min,
                time_resolution = device_capabilities.rydberg.global_value.time_resolution,
                max_time = device_capabilities.rydberg.global_value.time_max,
                max_value = device_capabilities.rydberg.global_value.detuning_max,
                min_value = device_capabilities.rydberg.global_value.detuning_min,
                max_slope = device_capabilities.rydberg.global_value.detuning_slew_rate_max,
                value_resolution = device_capabilities.rydberg.global_value.detuning_resolution
            ),
            δ = (
                min_time_step = device_capabilities.rydberg.local_value.time_delta_min,
                time_resolution = device_capabilities.rydberg.local_value.time_resolution,
                max_time = device_capabilities.rydberg.global_value.time_max,
                max_value = device_capabilities.rydberg.local_value.detuning_max,
                min_value = device_capabilities.rydberg.local_value.detuning_min,
                max_slope = device_capabilities.rydberg.local_value.detuning_slew_rate_max,
                value_resolution = device_capabilities.rydberg.local_value.common_detuning_resolution,
                local_mask_resolution = device_capabilities.rydberg.local_value.local_detuning_resolution
            )
        )
    end
end
