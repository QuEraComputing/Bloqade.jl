using Test
using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema:
    schema_parse_pwl_field,
    schema_parse_pwc_field,
    schema_parse_rydberg_fields


@testset "schema_parse_pwl_field" begin
    wf = piecewise_linear(;clocks=[0,2,3,4],values=[1,2,3,4])
    @test schema_parse_pwl_field(:wf,wf) == wf
    
    wf = Waveform(t->t^2,1)
    @test_throws ErrorException schema_parse_pwl_field(:wf,wf)
end

@testset "schema_parse_pwc_field" begin
    wf = piecewise_constant(;clocks=[0,2,3,4],values=[1,2,3])
    @test schema_parse_pwc_field(:wf,wf) == wf
    
    wf = Waveform(t->t^2,1)
    @test_throws ErrorException schema_parse_pwc_field(:wf,wf)
end

@testset "schema_parse_rydberg_fields" begin
    Ω = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,1,0])
    Δ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1,-1])
    ϕ = piecewise_constant(;clocks=Float64[0,1,2,3],values=Float64[0,1,-1])
    atoms = [(5.0*i,0.0) for i in 1:10]
    h = rydberg_h(atoms, Ω=Ω, Δ=Δ, ϕ=ϕ)

    @test schema_parse_rydberg_fields(h) == (atoms,ϕ,Ω,Δ,nothing,1.0)

    
    
end