
using BloqadeSchema
using BloqadeExpr
using Configurations
using Test
using BloqadeSchema
using BloqadeWaveforms
using Unitful: μs, s, MHz, rad 
    


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
    A = constant(;duration=2,value=1)
    for x in [0.1,1,1.5]
        B = constant(;duration=x,value=1)
        @test 2-x ≈ BloqadeSchema.norm_diff_durations(A,B)
    end
end


# test for all warning and error functions
@testset "warning and errors" begin
    for value in [nothing,1.0,rand(10),t->t^2]
        @test_throws ErrorException BloqadeSchema.check_waveform(value,:test)
    end
    
    wf = Waveform(t->t^2,π)
    wfs = [i*wf for i in rand(10)]

    @test_throws ErrorException BloqadeSchema.check_global(wfs,:test)

    time_res = 1.0e-3
    @test_logs (:info,"waveform test duration will be rounded during hardware transformation.") BloqadeSchema.warn_duration(time_res,wf,:test)
    @test_logs (:info,"waveform test duration will be rounded during hardware transformation.") BloqadeSchema.warn_duration(time_res,wfs,:test)

end

@testset "pin_waveform_edges" begin
    wf = Waveform(t->1+t^2,1)
    max_slope = 100
    begin_value = 0.5
    end_value = 1.0
    new_wf = BloqadeSchema.pin_waveform_edges(wf,max_slope,begin_value,end_value)
    @test new_wf(0.0) == begin_value
    @test new_wf(1.0) == end_value

    wf = Waveform(t->t^2,1)
    max_slope = 100
    begin_value = 0.0
    end_value = 0.0
    new_wf = BloqadeSchema.pin_waveform_edges(wf,max_slope,begin_value,end_value)
    @test new_wf(0.0) == begin_value
    @test new_wf(1.0) == end_value

    
    wf = Waveform(t->1-t^2,1)
    max_slope = 100
    begin_value = 0.0
    end_value = 0.0
    new_wf = BloqadeSchema.pin_waveform_edges(wf,max_slope,begin_value,end_value)
    @test new_wf(0.0) == begin_value
    @test new_wf(1.0) == end_value

end


@testset "clip_waveform" begin
    wf = piecewise_linear(;clocks=[0,1,2,3,4,5,6],values=[0,2,2,0,-2,-2,0])
    @test BloqadeSchema.clip_waveform(wf,:wf,-1,1) == piecewise_linear(;
        clocks=[0,1,2,3,4,5,6],
        values=[0,1,1,0,-1,-1,0])

    @test_logs (:info,"During hardware transform: wf(t) falls outside of hardware bounds, clipping to maximum/minimum.") BloqadeSchema.clip_waveform(wf,:wf,-1,1)
end

# @testset "find_local_masks" begin
#     Δ_values = 1 .- 2 .* rand(100,1)
#     Δ_mask = ones(1,10)
#     δ_values = rand(100,1)
#     δ_mask = rand(1,10)
#     δ_mask = (δ_mask .- minimum(δ_mask))/(maximum(δ_mask) .- minimum(δ_mask))

#     values = Δ_values * Δ_mask .+ δ_values * δ_mask

#     ((parsed_δ_values,parsed_δ_mask),(parsed_Δ_values,parsed_Δ_mask)) = BloqadeSchema.find_local_masks(values)

#     @test parsed_Δ_values ≈ reshape(Δ_values,(100,))
#     @test parsed_Δ_mask ≈ reshape(Δ_mask,(10,))
#     @test parsed_δ_values ≈ reshape(δ_values,(100,))
#     @test parsed_δ_mask ≈ reshape(δ_mask,(10,))
# end

# global test cases:
# 1. arbitrary waveform, check error. This includes boundary conditions for Ω.
# 2. piecewise linear waveform, check that the result doesn't change
@testset "hardware_transform_ϕ" begin
    device_capabilities = get_device_capabilities()

    wf = Waveform(t->1+t^2,1)
    wf_1,error_1 = BloqadeSchema.hardware_transform_Δ(wf,device_capabilities)
    wf_2,error_2 = BloqadeSchema.hardware_transform_Δ(wf_1,device_capabilities)

    @test wf_1 == wf_2
    @test error_2 == 0
end

@testset "hardware_transform_Ω" begin
    device_capabilities = get_device_capabilities()

    wf = Waveform(t->1+t^2,1)
    new_wf,error = BloqadeSchema.hardware_transform_Ω(wf,device_capabilities)
    @test new_wf(0.0) == 0
    @test new_wf(1.0) == 0

    wf = Waveform(t->t^2,1)
    new_wf,error = BloqadeSchema.hardware_transform_Ω(wf,device_capabilities)
    @test new_wf(0.0) == 0
    @test new_wf(1.0) == 0

    
    wf = Waveform(t->1-t^2,1)
    new_wf,error = BloqadeSchema.hardware_transform_Ω(wf,device_capabilities)
    @test new_wf(0.0) == 0
    @test new_wf(1.0) == 0
end

# local test case
# for piecewise linear waveforms Δ and δ, Δi that satisfy resolution 
# conditions check that this function can return a waveform that passes '≈' test. 
@testset "hardware_transform_Δ" begin
    device_capabilities = get_device_capabilities()

    wf = Waveform(t->1+t^2,1)
    wf_1,error_1 = BloqadeSchema.hardware_transform_Δ(wf,device_capabilities)
    wf_2,error_2 = BloqadeSchema.hardware_transform_Δ(wf_1,device_capabilities)

    @test wf_1 == wf_2
    @test error_2 == 0
end


@testset "hardware_transform_parse" begin

    atoms = [(1.14,4.0),(2.2,1),(5,3.37),(10,5.35)]
    dc = get_device_capabilities()
    Δ = piecewise_linear(;clocks=[0.0,1.0,2.0,3.0],values=[1.0,0.0,0.0,-1.0])
    Ω = piecewise_linear(;clocks=[0.0,1.0,2.0,3.0],values=[0.0,1.0,1.0,0.0])
    ϕ = piecewise_linear(;clocks=[0.0,1.0,2.0,3.0],values=[0.0,1,-1,0.0])
    H = rydberg_h(atoms,Δ=Δ,ϕ=ϕ,Ω=Ω)

    Δ_mask = (Δ=Δ,δ=nothing,Δi=1.0)

    (new_atoms,new_ϕ,new_Ω,new_Δ,info) = BloqadeSchema.hardware_transform_parse(H,dc)

    mse_atoms = (0.04+0.0+0.03+0.05)/4.0

    @test new_atoms == [(1.1,4.0),(2.2,1.0),(5.0,3.4),(10.0,5.3)]
    @test Δ == new_Δ
    @test Ω == new_Ω
    @test ϕ == new_ϕ
    @test info.Δ == 0
    @test info.Ω == 0
    @test info.ϕ == 0
    @test info.mse_atoms ≈ mse_atoms # up to rouding errors
    @test info.Δ_mask == Δ_mask

end
