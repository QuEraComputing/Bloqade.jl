using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema
using JSON



@testset "to_schema" begin
    T = 1
    s = T/10
    atoms = [(i,0) for i in 1:10]
    Δ = Waveform(t->cos(2π*t/T),T)
    # Δi=Δ
    amps = rand(length(atoms))
    Δi = [1.0*Δ for a in amps]

    Ω = Waveform(t->sin(2π*t/T)^2,T)
    ϕ = piecewise_linear(;clocks=[0,T/2-s,T/2+s,T],values=[0.0,0.0,π,π])

    H = rydberg_h(atoms;Ω=Ω,Δ=Δi,ϕ=ϕ)

    h = BloqadeSchema.to_json(H,waveform_tolerance=1e-2)
    H2 = BloqadeSchema.from_json(h)

    println(norm.(H.detuning_term.Δ .- H2.detuning_term.Δ))

end

