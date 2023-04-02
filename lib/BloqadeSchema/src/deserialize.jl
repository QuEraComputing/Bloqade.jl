function Configurations.from_dict(::Type{Lattice}, ::Type{NTuple{2,Float64}}, x)
    return (x[1], x[2])
end

function Configurations.from_dict(::Type{RabiFrequencyAmplitude}, d::Dict{String,<:Any})
    return RabiFrequencyAmplitude(Configurations.from_dict(GlobalField, d["global"]))
end

function Configurations.from_dict(::Type{RabiFrequencyPhase}, d::Dict{String,<:Any})
    return RabiFrequencyPhase(Configurations.from_dict(GlobalField, d["global"]))
end

function Configurations.from_dict(::Type{Detuning}, d::Dict{String,<:Any})
    return if haskey(d, "local")
        Detuning(
            global_value=Configurations.from_dict(GlobalField, d["global"]),
            local_value=Configurations.from_dict(LocalField, d["local"])
        )
    else
        Detuning(global_value = Configurations.from_dict(GlobalField, d["global"]), local_value = nothing)
    end
end

function Configurations.from_dict(::Type{RydbergCapabilities}, d::Dict{String,<:Any})
    return if haskey(d, "local")
        RydbergCapabilities(
            c6_coefficient=d["c6_coefficient"],
            global_value=Configurations.from_dict(RydbergGlobalCapabilities, d["global"]),
            local_value=Configurations.from_dict(RydbergLocalCapabilities, d["local"])
        )
    else
        RydbergCapabilities(c6_coefficient=d["c6_coefficient"],global_value = Configurations.from_dict(RydbergGlobalCapabilities, d["global"]), local_value = nothing)
    end
end
