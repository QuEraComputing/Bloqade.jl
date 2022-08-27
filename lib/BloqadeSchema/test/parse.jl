using Test
using BloqadeWaveforms
using BloqadeSchema



@testset "BloqadeSchema.parse_dynamic_rydberg_Ω" begin

    wf = Waveform(t->sin(2π*t)^2,2)
    Amps = [0.1*i for i in 1:10]
    wfs = [a*wf for a in Amps]
    duration = wf.duration
    # local drive test
    @test_throws ErrorException new_Amps,new_wf,duration = BloqadeSchema.parse_dynamic_rydberg_Ω(wfs) 
    # global drive test
    @test BloqadeSchema.parse_dynamic_rydberg_Ω(wf;duration=duration) == (1.0,wf,duration)

end

@testset "BloqadeSchema.parse_static_rydberg_Ω" begin
    duration = 4
    max_slope = 10.0
    min_step = 0.01
    Amps = [1,2,3,4]
    # local constant test
    @test_throws ErrorException new_Amps,new_wf = BloqadeSchema.parse_static_rydberg_Ω(Amps,duration,max_slope,min_step) 
    # global constant test
    wf = 1.0
    @test BloqadeSchema.parse_static_rydberg_Ω(wf,duration,max_slope,min_step) == (
        1.0,
        piecewise_linear(;clocks=Float64[0.0,0.1,duration-0.1,duration],values=Float64[0.0,wf,wf,0.0])
    )

end


@testset "BloqadeSchema.parse_dynamic_rydberg_Δ" begin

    # dynamic tests
    wf = Waveform(t->sin(2π*t)^2,2)
    Amps = [0.1*i for i in 1:10]
    wfs = [a*wf for a in Amps]

    # local drive test
    new_Amps,new_wf,duration = BloqadeSchema.parse_dynamic_rydberg_Δ(wfs) 
    new_values = sample_values(new_wf)

    @test duration == new_wf.duration
    for (i,wf_i) in enumerate(wfs)
        values = sample_values(wf_i)
        @test all(values .≈ new_Amps[i] .* new_values)
    end

    # global drive test
    @test BloqadeSchema.parse_dynamic_rydberg_Δ(wf;duration=duration) == (1.0,wf,duration)
end

@testset "BloqadeSchema.parse_static_rydberg_Δ" begin
    # constant tests
    duration = 4
    max_slope = 10.0
    min_step = 0.01
    Amps = [1,2,3,4]

    # local static test
    new_Amps,new_wf = BloqadeSchema.parse_static_rydberg_Δ(Amps,duration,max_slope,min_step)
    @test new_wf == piecewise_linear(;clocks=Float64[0,duration],values=Float64[1.0,1.0])
    @test all(new_Amps .== Amps)        

    # constant value test
    wf = 1.0
    @test BloqadeSchema.parse_static_rydberg_Δ(wf,duration,max_slope,min_step) == (
        1.0,
        piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[wf,wf])
    )
end
