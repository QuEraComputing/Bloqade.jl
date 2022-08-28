using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema
using JSON



@testset "to_schema" begin
    T = 2.0
    atoms = [(i,0) for i in 1:10]
    Δ = Waveform(t->cos(2π*t/T),T)
    # Δi=Δ
    Amps = rand(length(atoms))
    Δi = [a*Δ for a in Amps]

    Ω = Waveform(t->sin(2π*t/T)^2,T)
    ϕ = piecewise_constant(;clocks=[0,T/2,T],values=[0,π])
    # ϕ = Waveform(t->π*sin(π*t/2/T)^2,T)

    H = rydberg_h(atoms;Ω=Ω,Δ=Δi,ϕ=ϕ)

    h = BloqadeSchema.to_json(H,waveform_tolerance=1e-1)
    H2 = BloqadeSchema.from_json(h)

    println(H2)
end

