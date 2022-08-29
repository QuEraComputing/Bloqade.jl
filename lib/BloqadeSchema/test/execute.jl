using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema
using JSON



@testset "to_schema" begin
    T = 2.0
    s = 0.2
    atoms = [(i,0) for i in 1:10]
    Δ = Waveform(t->cos(2π*t/T),T)
    # Δi=Δ
    Amps = rand(length(atoms))
    Δi = [a*Δ for a in Amps]

    Ω = Waveform(t->sin(2π*t/T)^2,T)
    ϕ = piecewise_linear(;clocks=[0,T/2-s,T/2+s,T],values=[0.0,0.0,π,π])

    H = rydberg_h(atoms;Ω=Ω,Δ=Δi,ϕ=ϕ)

    h = BloqadeSchema.to_json(H,waveform_tolerance=1e-3)
    H2 = BloqadeSchema.from_json(h)
    
end

