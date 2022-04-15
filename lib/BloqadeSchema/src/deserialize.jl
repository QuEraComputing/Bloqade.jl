function Configurations.from_dict(::Type{Lattice}, ::Type{NTuple{2, Float64}}, x)
    return (x[1], x[2])
end

function Configurations.from_dict(::Type{RydbergRabiFrequencyAmplitude}, d::Dict{String, Any})
    return RydbergRabiFrequencyAmplitude(from_dict(RydbergRabiFrequencyAmplitudeGlobal, d["global"]))
end

function Configurations.from_dict(::Type{RydbergRabiFrequencyPhase}, d::Dict{String, Any})
    return RydbergRabiFrequencyPhase(from_dict(RydbergRabiFrequencyPhaseGlobal, d["global"]))
end

function Configurations.from_dict(::Type{RydbergDetuning}, d::Dict{String, Any})
    return if haskey(d, "local")
        RydbergDetuning(
            from_dict(RydbergDetuningGlobal, d["global"]),
            from_dict.(RydbergDetuningLocal, d["local"]),
        )
    else
        RydbergDetuning(
            from_dict(RydbergDetuningGlobal, d["global"]),
            nothing,
        )
    end
end


# BloqadeSchema.RydbergHamiltonian
# from_dict(::Type{OptionType}, d::AbstractDict{String}; kw...) where {OptionType}