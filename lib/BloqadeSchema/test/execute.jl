using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema
using JSON


@testset "to_schema" begin
    T = 4
    atoms = [(i,0) for i in 1:10]

    values = [1.0,nothing,Waveform(t->sin(2π*t/T)^2,T/2)]
    for Ω in values
        for ϕ in values
            for Δ in values
                if any((f isa BloqadeSchema.DynamicParam) for f in [ϕ,Ω,Δ])
                    H = rydberg_h(atoms;Ω=Ω,Δ=Δ,ϕ=ϕ)
                    h = BloqadeSchema.to_json(H,waveform_tolerance=1e-2,warn=true)
                end
            end
        end
    end

end

