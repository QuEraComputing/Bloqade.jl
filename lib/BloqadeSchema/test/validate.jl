
using BloqadeWaveforms
using BloqadeSchema
using BloqadeExpr

@testset "run" begin
    atoms = [(i,0) for i in 1:10]
    Δ = Waveform(t->cos(t),1)
    Ω = Waveform(t->sin(t),1)
    ϕ = constant(;value=0,duration=1)
    h = rydberg_h(atoms,ϕ=ϕ,Ω=Ω,Δ=Δ)


    h_hardware,info = hardware_transform(h) 


    validate(h_hardware;warn=true)
end

