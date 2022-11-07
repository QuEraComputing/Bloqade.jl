function Configurations.to_dict(
    ::Type{T},
    x::Union{RydbergRabiFrequencyAmplitude,RydbergRabiFrequencyPhase},
    option::Configurations.ToDictOption,
) where {T}
    return Dict("global" => Configurations.to_dict(x.global_value))
end

function Configurations.to_dict(::Type{T}, x::Union{RydbergDetuning,RydbergCapabilities}, option::Configurations.ToDictOption) where {T}
    if isnothing(x.local_value)
        return Dict("global" => Configurations.to_dict(x.global_value))
    end

    return Dict("global" => Configurations.to_dict(x.global_value), "local" => Configurations.to_dict(x.local_value))
end


function Configurations.to_dict(::Type{T}, x::Lattice, option::Configurations.ToDictOption) where {T}
    return Dict("sites" => map(collect, x.sites), "filling" => x.filling)
end
