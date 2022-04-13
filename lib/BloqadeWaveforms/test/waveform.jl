using Test
using Intervals
using BloqadeWaveforms

@testset "Waveform" begin
    waveform = Waveform(t->2.2sin(t), duration=4π)
    @test waveform(0.1) ≈ 2.2 * sin(0.1)
    @test_throws ArgumentError waveform(0.1+4π)

    # constant bindings
    wf = sinusoidal(;duration=4π, amplitude=2.2)
    @test wf(0.1) ≈ 2.2 * sin(0.1)

    show(stdout, wf)
    println(stdout)
    show(stdout, MIME"text/plain"(), wf)
    println(stdout)

    wf = constant(;duration=4π, value=2.2)
    @test wf(0.1) ≈ 2.2

    show(stdout, wf)
    println(stdout)
    show(stdout, MIME"text/plain"(), wf)
    println(stdout)
end

@testset "linear_ramp" begin
    waveform = linear_ramp(duration=0.5, start_value=0.0, stop_value=1.0)
    @test waveform(0.1) ≈ 2 * 0.1
    @test_throws ArgumentError waveform(0.8)
end

@testset "constant" begin
    waveform = constant(duration=0.5, value=2.1)
    @test waveform(0.1) ≈ 2.1
    @test_throws ArgumentError waveform(0.6)
end

@testset "piecewise" begin
    @testset "piecewise_constant" begin
        waveform = piecewise_constant(clocks=[0.0, 0.2, 0.5], values=[0.0, 1.5, 3.1])
        @test waveform(0.0) ≈ 0.0
        @test waveform(0.1) ≈ 0.0
        @test waveform(0.2) ≈ 1.5
        @test waveform(0.3) ≈ 1.5
        @test waveform(0.5) ≈ 3.1
        @test_throws ArgumentError waveform(0.6) ≈ 3.1
    
        waveform = piecewise_constant(clocks=[0.0, 0.2, 0.5], values=[0.0, 1.5, 3.1], duration=1.1)
        @test waveform(0.6) ≈ 3.1
    end
    
    @testset "piecewise_linear" begin
        waveform = piecewise_linear(clocks=[0.0, 0.2, 0.5, 0.8, 1.0], values=[0.0, 1.5, 3.1, 3.1, 0.0])
        @test waveform(0.1) ≈ 0.75
        @test waveform(0.6) ≈ 3.1
        @test waveform(1.0) ≈ 0.0
        @test_throws ArgumentError waveform(1.1)
    end
end

@testset "waveform + waveform" begin
    wf1 = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0)
    wf2 = Waveform(sin, duration=2.2)
    wf3 = wf1 + wf2
    @test wf3(0.1) ≈ wf1(0.1) + wf2(0.1)

    # sum + other
    wfp = piecewise_constant(clocks=[0.0, 0.3, 0.5], values=[0.0, 1.1, 0.5], duration=2.2)
    wf4 = wf3 + wfp
    wf5 = wfp + wf3
    @test wf4 isa Waveform
    @test wf5 isa Waveform

    # sum + sum
    wf5 = wf3 + wf3
    @test wf5 isa Waveform

    wf1 = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0)
    wf2 = Waveform(sin, duration=2.1)
    @test_throws ArgumentError wf1 + wf2
end

@testset "-waveform" begin
    wf1 = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0)
    wf2 = Waveform(sin, duration=2.2)
    wf3 = -wf1
    @test wf3(0.1) ≈ -wf1(0.1)

    wf4 = wf1 - wf2
    @test wf4(0.1) ≈ wf1(0.1) - wf2(0.1)

    wf5 = -wf1 + wf2
    @test wf5(0.1) ≈ -wf1(0.1) + wf2(0.1)

    wf6 = -wf1 - wf2
    @test wf6(0.1) ≈ -wf1(0.1) - wf2(0.1)

    wf7 = wf2 - (-wf1)
    @test wf7(0.1) ≈ wf2(0.1) + wf1(0.1)
end

@testset "alpha*waveform" begin
    wf1 = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0)
    wf2 = 2.0 * wf1
    @test wf2(0.1) ≈ 2 * wf1(0.1)

    wf3 = wf1 * 2.0
    @test wf3(0.1) ≈ wf1(0.1) * 2
end

@testset "alpha*waveform" begin
    wf1 = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0)
    wf2 = 2.0 / wf1
    @test wf2(0.1) ≈ 2 / wf1(0.1)

    wf3 = wf1 / 2.0
    @test wf3(0.1) ≈ wf1(0.1) / 2
end

@testset "broadcast" begin
    wf1 = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0)
    wf2, wf3 = [2.0, 3.0] .* wf1
    @test wf2(0.1) ≈ 2.0 * wf1(0.1)
    @test wf3(0.1) ≈ 3.0 * wf1(0.1)
end

@testset "append(waveforms...)" begin
    wf1 = Waveform(sin, duration=2.2)
    wf2 = linear_ramp(;start_value=0.0, stop_value=1.1, duration=0.5)
    waveform = append(wf1, wf2)

    @testset "(::CompositeWaveform)(t::Real)" begin
        @test waveform(0.1) ≈ wf1(0.1)
        @test waveform(2.5) ≈ wf2(0.3)
        @test_throws ArgumentError waveform(5.2)
    end

    @testset "sample_values(::CompositeWaveform, $dt)" for dt in [1e-3, 1e-3, 1e-4]
        values = sample_values(waveform; dt)
        @test length(values) == length(sample_clock(waveform; dt))
        wf_value_1 = wf1.(sample_clock(wf1; dt))
        wf_value_2 = wf2.(sample_clock(wf2; dt))

        @test values[1:length(wf_value_1)] ≈ wf_value_1[1:end]
        @test values[length(wf_value_1)+1:end] ≈ wf_value_2[2:end]
    end
end

@testset "slicing" begin
    wf = Waveform(sin, duration=4π)
    wfs = wf[0.5..2π]
    @test wfs.duration == 2π-0.5
    @test wfs(0.0) ≈ wf(0.5)
    @test wfs(0.1) ≈ wf(0.6)

    @test_throws ArgumentError wf[0.5..(4π+2)]
end

@testset "piecewise_constant/linear assertions" begin
    @test_throws ArgumentError piecewise_constant(clocks=[0,1], values=[1,2,3])
    @test_throws ArgumentError piecewise_constant(clocks=[-1,1], values=[1,2,3])
    @test_throws ArgumentError piecewise_linear(clocks=[0,1], values=[1,2,3])
    @test_throws ArgumentError piecewise_linear(clocks=[-1,1], values=[1,2,3])
end
