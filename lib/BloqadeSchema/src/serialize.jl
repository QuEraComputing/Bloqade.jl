function Configurations.to_dict(::Type{T}, x::RydbergRabiFrequencyAmplitude, option::Configurations.ToDictOption) where T
    return Dict(
        "global" => to_dict(x.global_value)
    )
end
