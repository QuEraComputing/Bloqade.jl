using Test
using BloqadeWaveforms
using BloqadeSchema: 
    parse_dynamic_rydberg_Ω,
    parse_dynamic_rydberg_Δ,
    parse_dynamic_rydberg_ϕ,
    parse_static_rydberg_Ω,
    parse_static_rydberg_Δ,
    parse_static_rydberg_ϕ,
    convert_units,
    TaskSpecification

using Unitful: μs, s, MHz, rad 
    

@testset "convert_units" begin
    T = 1 # seconds
    T_list = [i for i in 1:10]
    pair = (1,2)
    pair_list = [(i,j) for i in 1:5 for j in 1:5]
    @test 1e6 ≈ convert_units(T,s,μs)
    @test (1e6 .* T_list) ≈ convert_units(T_list,s,μs)
    @test all((1e6 .* pair ).≈ convert_units.(pair,s,μs))
    # @test all([1e6.*p for p in pair_list] .≈ convert_units.(pair_list,s,μs))

end

@testset "parse_dynamic_rydberg_Ω" begin



end

@testset "parse_static_rydberg_Ω" begin
  

end

@testset "parse_dynamic_rydberg_ϕ" begin

end

@testset "parse_static_rydberg_ϕ" begin
   

end


@testset "parse_dynamic_rydberg_Δ" begin

end

@testset "parse_static_rydberg_Δ" begin

end
