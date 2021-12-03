using Test
using Intervals
using EaRydWaveforms
using EaRydWaveforms: SumWaveform, NegWaveform, SlicedWaveform

@testset "FunctionWaveform" begin
    waveform = FunctionWaveform(t->2.2sin(t), duration=4π)
    @test waveform(0.1) ≈ 2.2 * sin(0.1)
    @test_throws ArgumentError waveform(0.1+4π)

    # constant bindings
    wf = SinusoidalWaveform(;duration=4π, amplitude=2.2)
    @test wf(0.1) ≈ 2.2 * sin(0.1)

    show(stdout, wf)
    println(stdout)
    show(stdout, MIME"text/plain"(), wf)
    println(stdout)

    wf = ConstantWaveform(;duration=4π, amplitude=2.2)
    @test wf(0.1) ≈ 2.2

    show(stdout, wf)
    println(stdout)
    show(stdout, MIME"text/plain"(), wf)
    println(stdout)
end

@testset "RampWaveform" begin
    waveform = RampWaveform(duration=0.5, start=0.0, stop=1.0)
    @test waveform(0.1) ≈ 2 * 0.1
    @test_throws ArgumentError waveform(0.8)
end

@testset "ConstantWaveform" begin
    waveform = ConstantWaveform(duration=0.5, amplitude=2.1)
    @test waveform(0.1) ≈ 2.1
    @test_throws ArgumentError waveform(0.6)
end

@testset "Piecewise" begin
    @testset "PiecewiseConsant" begin
        waveform = PiecewiseConsant(clocks=[0.0, 0.2, 0.5], values=[0.0, 1.5, 3.1])
        @test waveform(0.0) ≈ 0.0
        @test waveform(0.1) ≈ 0.0
        @test waveform(0.2) ≈ 1.5
        @test waveform(0.3) ≈ 1.5
        @test waveform(0.5) ≈ 3.1
        @test_throws ArgumentError waveform(0.6) ≈ 3.1
    
        waveform = PiecewiseConsant(clocks=[0.0, 0.2, 0.5], values=[0.0, 1.5, 3.1], duration=1.1)
        @test waveform(0.6) ≈ 3.1
    end
    
    @testset "PiecewiseLinear" begin
        waveform = PiecewiseLinear(clocks=[0.0, 0.2, 0.5, 0.8, 1.0], values=[0.0, 1.5, 3.1, 3.1, 0.0])
        @test waveform(0.1) ≈ 0.75
        @test waveform(0.6) ≈ 3.1
        @test waveform(1.0) ≈ 0.0
        @test_throws ArgumentError waveform(1.1)
    end
end

@testset "SumWaveform" begin
    wf1 = RampWaveform(;duration=2.2, start=0.0, stop=1.0)
    wf2 = FunctionWaveform(sin, duration=2.2)
    wf3 = wf1 + wf2
    @test wf3(0.1) ≈ wf1(0.1) + wf2(0.1)

    # sum + other
    wfp = PiecewiseConsant(clocks=[0.0, 0.3, 0.5], values=[0.0, 1.1, 0.5], duration=2.2)
    wf4 = wf3 + wfp
    wf5 = wfp + wf3
    @test wf4 isa SumWaveform
    @test wf5 isa SumWaveform
    @test length(wf4.waveforms) == 3
    @test length(wf5.waveforms) == 3
    @test wf4.waveforms isa Tuple{RampWaveform, FunctionWaveform, PiecewiseConsant}
    @test wf5.waveforms isa Tuple{PiecewiseConsant, RampWaveform, FunctionWaveform}

    # sum + sum
    wf5 = wf3 + wf3
    @test wf5 isa SumWaveform
    @test length(wf5.waveforms) == 4

    wf1 = RampWaveform(;duration=2.2, start=0.0, stop=1.0)
    wf2 = FunctionWaveform(sin, duration=2.1)
    @test_throws ArgumentError wf1 + wf2
end

@testset "NegWaveform" begin
    wf1 = RampWaveform(;duration=2.2, start=0.0, stop=1.0)
    wf2 = FunctionWaveform(sin, duration=2.2)
    wf3 = -wf1
    @test wf3 isa NegWaveform
    @test wf3(0.1) ≈ -wf1(0.1)
    @test -wf3 isa typeof(wf1)

    wf4 = wf1 - wf2
    @test wf4 isa SumWaveform
    @test wf4(0.1) ≈ wf1(0.1) - wf2(0.1)

    wf5 = -wf1 + wf2
    @test wf5(0.1) ≈ -wf1(0.1) + wf2(0.1)

    wf6 = -wf1 - wf2
    @test wf6(0.1) ≈ -wf1(0.1) - wf2(0.1)

    wf7 = wf2 - (-wf1)
    @test wf7(0.1) ≈ wf2(0.1) + wf1(0.1)
end

@testset "CompositeWaveform" begin
    waveform = CompositeWaveform(
        FunctionWaveform(sin, duration=2.2),
        RampWaveform(;duration=0.5, start=0.0, stop=1.0)
    )

    @testset "(::CompositeWaveform)(t::Real)" begin
        @test waveform(0.1) ≈ waveform[1](0.1)
        @test waveform(2.5) ≈ waveform[2](0.3)
        @test_throws ArgumentError waveform(5.2)
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

@testset "SlicedWaveform" begin
    wf = FunctionWaveform(sin, duration=4π)
    wfs = wf[0.5..2π]
    @test duration(wfs) ≈ span(0.5..2π)
    @test wfs(0.0) ≈ wf(0.5)
    @test wfs(0.1) ≈ wf(0.6)

    @test_throws ArgumentError wf[0.5..(4π+2)]
end
