using Test
# using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema:
    get_rydberg_params,
    schema_parse_field

using BloqadeWaveforms


@testset "get_rydberg_params" begin
    atoms = [(i,i) for i in 1:10]

    values = [nothing,1,t->t^2,rand(10),[t->rand()*t for i in 1:10]]

    for ϕ in values, Ω in values, Δ in values
        h = rydberg_h(atoms,ϕ=ϕ,Ω=Ω,Δ=Δ)
        # catches the this weird edge case 
        Ω = (isnothing(Ω) && !isnothing(ϕ) ? 0 : Ω)
        @test (atoms,ϕ,Ω,Δ) == BloqadeSchema.get_rydberg_params(h)
    end
end

@testset "schema_parse_field" begin
    wf = piecewise_linear(;clocks=[0,2,3,4],values=[1,2,3,4])
    @test schema_parse_field(:wf,wf) == wf
    
    wf = Waveform(t->t^2,1)
    @test_throws ErrorException schema_parse_field(:wf,wf)
end

"""
@testset "schema_parse_ϕ" begin
    wf = piecewise_linear(;clocks=[0,2,3,4],values=[1,2,3,4])
    @test schema_parse_ϕ(wf) == wf
    @test_throws ErrorException schema_parse_ϕ(Waveform(t->t^2,1))
end

@testset "schema_parse_Ω" begin
    wf = piecewise_linear(;clocks=[0,2,3,4],values=[1,2,3,4])
    @test schema_parse_Ω(wf) == wf
    @test_throws ErrorException schema_parse_Ω(Waveform(t->t^2,1))
end

@testset "schema_parse_Δ" begin
    clocks = Float64[0,1,2,3]
    Δ_values = Float64[1,0,0,-1]
    δ_values = Float64[1,0,0,1]
    Δi = rand(10,)
    Δi = (Δi .- minimum(Δi))./(maximum(Δi)-minimum(Δi))
    Δi = set_resolution.(Δi,0.001)
    Δ = piecewise_linear(;clocks=clocks,values=Δ_values)
    δ = piecewise_linear(;clocks=clocks,values=δ_values)
    @test (Δ,nothing,1.0) == schema_parse_Δ(Δ)
    @test_throws ErrorException schema_parse_Δ(Waveform(t->t^2,1.0))
    @test_throws ErrorException schema_parse_Δ([Waveform(t->t^2,1.0) for i in 1:4])

    local_Δ = [piecewise_linear(;clocks=clocks,values=Δ_values+δi.*δ_values) for δi in Δi]

    # @test (Δ,δ,Δi) == schema_parse_Δ(local_Δ)
    (parsed_Δ,parsed_δ,parsed_Δi) = schema_parse_Δ(local_Δ)
    println(Δ.f.values.-parsed_Δ.f.values)
    @test parsed_Δ == Δ
    @test parsed_δ == δ
    @test parsed_Δi == Δi
end
"""