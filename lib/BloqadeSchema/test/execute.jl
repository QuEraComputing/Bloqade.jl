using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema
using JSON



@testset "to_schema" begin
    T = 1
    atoms = [(i,0) for i in 1:10]
    Δ = Waveform(t->cos(2π*t/T),T)
    Amps = [0.1*i for i in 1:10]
    Δi = [a*Δ for a in Amps]

    Ω = Waveform(t->sin(2π*t/T)^2,T)

    H = rydberg_h(atoms;Ω=Ω,Δ=Δi)

    h = BloqadeSchema.to_json(H,waveform_tolerance=1e-2)
    JSON.print(JSON.parse(h), 4)
end

