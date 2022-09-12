using Test
# using BloqadeWaveforms
using BloqadeExpr
using BloqadeSchema:
    get_rydberg_params,
    schema_parse_ϕ,
    schema_parse_Ω,
    schema_parse_Δ

using BloqadeWaveforms


@testset "schema_parsing" begin
    atoms = [(i,i) for i in 1:10]

    values = [1,t->t^2,rand(10),[t->rand()*t for i in 1:10]]

    for ϕ in values, Ω in values, Δ in values
        h = rydberg_h(atoms,ϕ=ϕ,Ω=Ω,Δ=Δ)
        @test (atoms,ϕ,Ω,Δ) == BloqadeSchema.get_rydberg_params(h)

    end

    wf = piecewise_linear(;clocks=[0,2,3,4],values=[1,2,3,4])
    @test schema_parse_ϕ(wf) == wf
    @test schema_parse_Ω(wf) == wf
    @test schema_parse_Δ(wf) == (wf,nothing,1.0)

end
