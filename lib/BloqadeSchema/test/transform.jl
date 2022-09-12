

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

@testset "norm_diff_durations" begin
    
end


# test for all warning and error functions
@testset "warning and errors" begin

end

@testset "pin_waveform_edges" begin
    wf = Waveform(t->1+t^2,1)
    max_slope = 100
    begin_value = 0.5
    end_value = 1.0
    new_wf = BloqadeSchema.pin_waveform_edges(wf,max_slope,begin_value,end_value)
    @test new_wf(0.0) == begin_value
    @test new_wf(1.0) == end_value

end

@testset "find_local_masks" begin
    Δ_values = 1 .- 2 .* rand(100,1)
    Δ_mask = ones(1,10)
    δ_values = rand(100,1)
    δ_mask = rand(1,10)
    δ_mask = (δ_mask .- minimum(δ_mask))/(maximum(δ_mask) .- minimum(δ_mask))

    values = Δ_values * Δ_mask .+ δ_values * δ_mask

    ((parsed_Δ_values,parsed_Δ_mask),(parsed_δ_values,parsed_δ_mask)) = BloqadeSchema.find_local_masks(values)
    @test parsed_Δ_values ≈ Δ_values
    @test parsed_Δ_mask ≈ Δ_mask
    @test parsed_δ_values ≈ δ_values
    @test parsed_δ_mask ≈ δ_mask
end

# global test cases:
# 1. arbitrary waveform, check error. This includes boundary conditions for Ω.
# 2. piecewise linear waveform, check that the result doesn't change
@testset "hardware_transform_ϕ" begin
    
end

@testset "hardware_transform_Ω" begin
    
end

# local test case
# for piecewise linear waveforms Δ and δ, Δi that satisfy resolution 
# conditions check that this function can return a waveform that passes '≈' test. 
@testset "hardware_transform_Δ" begin
    
end