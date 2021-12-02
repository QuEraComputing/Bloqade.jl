using Test
using EaRydWaveforms: sample_values

clocks = collect(0:0.1:1)
values = cos.(clocks)

waveform = InterpolatedWaveform(clocks, values)
waveform = SinusoidalWaveform(2.2)

@testset "RampWaveform" begin
    waveform = RampWaveform(0.5, 0.0, 1.0)
     
end

@testset "CompositeWaveform" begin
    waveform = CompositeWaveform(
        SinusoidalWaveform(2.2),
        RampWaveform(0.5, 0.0, 1.0)
    )

    @testset "(::CompositeWaveform)(t::Real)" begin
        @test waveform(0.1) ≈ waveform[1](0.1)
        @test waveform(2.5) ≈ waveform[2](0.3)
    end

    @testset "sample_values(::CompositeWaveform, $dt)" for dt in [1e-3, 1e-3, 1e-4]
        values = sample_values(waveform, dt)
        @test length(values) == length(0:dt:duration(waveform))
        wf_value_1 = waveform[1].(0:dt:duration(waveform[1]))
        wf_value_2 = waveform[2].(0:dt:duration(waveform[2]))

        @test values[1:length(wf_value_1)-1] ≈ wf_value_1[1:end-1]
        @test values[length(wf_value_1):end] ≈ wf_value_2
    end        
end
