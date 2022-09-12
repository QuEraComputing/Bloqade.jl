

using BloqadeSchema
using Unitful: μs, s, MHz, rad 
    

@testset "convert_units" begin
    T = 1 # seconds
    T_list = [i for i in 1:10]
    pair = (1,2)
    pair_list = [(i,j) for i in 1:5 for j in 1:5]
    @test 1e6 ≈ BloqadeSchema.convert_units(T,s,μs)
    @test (1e6 .* T_list) ≈ BloqadeSchema.convert_units(T_list,s,μs)
    @test all((1e6 .* pair ).≈ BloqadeSchema.convert_units.(pair,s,μs))
    # @test all([1e6.*p for p in pair_list] .≈ convert_units.(pair_list,s,μs))

end

@testset "set_resolution" begin

    examples = [
        (1,10.123,10.0),
        (0.1,10.1256,10.1),
        (0.01,10.1256,10.13),
    ]
    for (δ,val,res) in examples
        @test BloqadeSchema.set_resolution(val,δ) == res
    end

end