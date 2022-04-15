function Configurations.to_dict(::Type{T}, x::Union{RydbergRabiFrequencyAmplitude,RydbergRabiFrequencyPhase}, option::Configurations.ToDictOption) where {T}
    return Dict(
        "global" => to_dict(x.global_value),
    )
end

function Configurations.to_dict(::Type{T}, x::RydbergDetuning, option::Configurations.ToDictOption) where {T}
    return Dict(
        "global" => to_dict(x.global_value),
        "local" => x.local_value
    )
end
