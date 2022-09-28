function Configurations.from_dict(::Type{Lattice}, ::Type{NTuple{2,Float64}}, x)
    return (x[1], x[2])
end

function Configurations.from_dict(::Type{RydbergRabiFrequencyAmplitude}, d::Dict{String,<:Any})
    return RydbergRabiFrequencyAmplitude(Configurations.from_dict(RydbergRabiFrequencyAmplitudeGlobal, d["global"]))
end

function Configurations.from_dict(::Type{RydbergRabiFrequencyPhase}, d::Dict{String,<:Any})
    return RydbergRabiFrequencyPhase(Configurations.from_dict(RydbergRabiFrequencyPhaseGlobal, d["global"]))
end

function Configurations.from_dict(::Type{RydbergDetuning}, d::Dict{String,<:Any})
    return if haskey(d, "local")
        RydbergDetuning(
            global_value=Configurations.from_dict(RydbergDetuningGlobal, d["global"]),
            local_value=Configurations.from_dict(RydbergDetuningLocal, d["local"])
        )
    else
        RydbergDetuning(global_value = Configurations.from_dict(RydbergDetuningGlobal, d["global"]), local_value = nothing)
    end
end
