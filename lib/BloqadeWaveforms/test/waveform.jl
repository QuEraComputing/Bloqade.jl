using Test
using Unitful: ms, rad
using Intervals
using BloqadeWaveforms

@testset "Waveform" begin
    waveform = Waveform(t -> 2.2sin(t), duration = 4π)
    @test waveform(0.1) ≈ 2.2 * sin(0.1)
    @test_throws ArgumentError waveform(0.1 + 4π)

    # constant bindings
    wf = sinusoidal(; duration = 2, amplitude = 2.2)
    @test wf(0.1) ≈ 2.2 * sin(2π * 0.1)

    show(stdout, wf)
    println(stdout)
    show(stdout, MIME"text/plain"(), wf)
    println(stdout)

    wf = constant(; duration = 4π, value = 2.2)
    @test wf(0.1) ≈ 2.2

    show(stdout, wf)
    println(stdout)
    show(stdout, MIME"text/plain"(), wf)
    println(stdout)
end

@testset "linear_ramp" begin
    waveform = linear_ramp(duration = 0.5, start_value = 0.0, stop_value = 1.0)
    @test waveform(0.1) ≈ 2 * 0.1
    @test_throws ArgumentError waveform(0.8)
end

@testset "constant" begin
    waveform = constant(duration = 0.5, value = 2.1)
    @test waveform(0.1) ≈ 2.1
    @test_throws ArgumentError waveform(0.6)
end

@testset "piecewise" begin
    @testset "piecewise_constant" begin
        waveform = piecewise_constant(clocks = [0.0, 0.2, 0.4, 0.5], values = [0.0, 1.5, 3.1])
        @test waveform(0.0) ≈ 0.0
        @test waveform(0.1) ≈ 0.0
        @test waveform(0.2) ≈ 1.5
        @test waveform(0.3) ≈ 1.5
        @test waveform(0.5) ≈ 3.1
        @test_throws ArgumentError waveform(0.6) ≈ 3.1

        waveform = piecewise_constant(clocks = [0.0, 0.2, 0.5, 1.1], values = [0.0, 1.5, 3.1])
        @test waveform(0.6) ≈ 3.1

        wf = piecewise_constant(clocks = [0.0ms, 0.1ms, 0.2ms, 0.3ms], values = [0.1, 1.1, 2.1] .* (rad / ms))
        @test wf(0.0) ≈ 0.0001
        @test wf(100.0) ≈ 0.0011
        @test wf(200.0) ≈ 0.0021
        @test wf.duration ≈ 300.0
    end

    @testset "piecewise_linear" begin
        waveform = piecewise_linear(clocks = [0.0, 0.2, 0.5, 0.8, 1.0], values = [0.0, 1.5, 3.1, 3.1, 0.0])
        @test waveform(0.1) ≈ 0.75
        @test waveform(0.6) ≈ 3.1
        @test waveform(1.0) ≈ 0.0
        @test_throws ArgumentError waveform(1.1)

        wf = piecewise_linear(clocks = [0.0ms, 0.1ms, 0.2ms], values = [0.1, 1.1, 2.1] .* (rad / ms))
        @test wf(0.0) ≈ 0.0001
        @test wf(100.0) ≈ 0.0011
        @test wf(200.0) ≈ 0.0021
    end

    @testset "promote_type" begin
        wf = piecewise_linear(clocks=[0,1.23], values= [0, 0])
        @test eltype(wf) === Float64
        wf = piecewise_constant(clocks=[0,1.23,1.5], values= [0, 0])
        @test eltype(wf) === Float64
    end
end

@testset "waveform + waveform" begin
    wf1 = linear_ramp(; duration = 2.2, start_value = 0.0, stop_value = 1.0)
    wf2 = Waveform(sin, duration = 2.2)
    wf3 = wf1 + wf2
    @test wf3(0.1) ≈ wf1(0.1) + wf2(0.1)

    # sum + other
    wfp = piecewise_constant(clocks = [0.0, 0.3, 0.5, 2.2], values = [0.0, 1.1, 0.5])
    wf4 = wf3 + wfp
    wf5 = wfp + wf3
    @test wf4 isa Waveform
    @test wf5 isa Waveform

    # sum + sum
    wf5 = wf3 + wf3
    @test wf5 isa Waveform

    wf1 = linear_ramp(; duration = 2.2, start_value = 0.0, stop_value = 1.0)
    wf2 = Waveform(sin, duration = 2.1)
    @test_throws ArgumentError wf1 + wf2
end

@testset "-waveform" begin
    wf1 = linear_ramp(; duration = 2.2, start_value = 0.0, stop_value = 1.0)
    wf2 = Waveform(sin, duration = 2.2)
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
    wf1 = linear_ramp(; duration = 2.2, start_value = 0.0, stop_value = 1.0)
    wf2 = 2.0 * wf1
    @test wf2(0.1) ≈ 2 * wf1(0.1)

    wf3 = wf1 * 2.0
    @test wf3(0.1) ≈ wf1(0.1) * 2
end

@testset "alpha*waveform" begin
    wf1 = linear_ramp(; duration = 2.2, start_value = 0.0, stop_value = 1.0)
    wf2 = 2.0 / wf1
    @test wf2(0.1) ≈ 2 / wf1(0.1)

    wf3 = wf1 / 2.0
    @test wf3(0.1) ≈ wf1(0.1) / 2
end

@testset "broadcast" begin
    wf1 = linear_ramp(; duration = 2.2, start_value = 0.0, stop_value = 1.0)
    wf2, wf3 = [2.0, 3.0] .* wf1
    @test wf2(0.1) ≈ 2.0 * wf1(0.1)
    @test wf3(0.1) ≈ 3.0 * wf1(0.1)
end

@testset "append(waveforms...)" begin
    wf1 = Waveform(sin, duration = 2.2)
    wf2 = linear_ramp(; start_value = 0.0, stop_value = 1.1, duration = 0.5)
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
    wf = Waveform(sin, duration = 4π)
    wfs = wf[0.5..2π]
    @test wfs.duration == 2π - 0.5
    @test wfs(0.0) ≈ wf(0.5)
    @test wfs(0.1) ≈ wf(0.6)

    @test_throws ArgumentError wf[0.5..(4π+2)]
end

@testset "piecewise_constant/linear assertions" begin
    @test_throws ArgumentError piecewise_constant(clocks = [0, 1], values = [1, 2, 3])
    @test_throws ArgumentError piecewise_constant(clocks = [-1, 1], values = [1, 2, 3])
    @test_throws ArgumentError piecewise_linear(clocks = [0, 1], values = [1, 2, 3])
    @test_throws ArgumentError piecewise_linear(clocks = [-1, 1], values = [1, 2, 3])
end

@testset "piecewise_constant/linear equality" begin
    @test piecewise_linear(clocks=[0, 3, 5], values=[3, 4, 5]) == piecewise_linear(clocks=[0.0, 3.0, 5], values=[3, 4, 5.0])
    @test piecewise_linear(clocks=[0, 3, 5], values=[3, 4, 8]) != piecewise_linear(clocks=[0.0, 3.0, 5], values=[3, 4, 5.0])

    @test piecewise_constant(clocks=[0, 3, 5, 6], values=[3, 4, 5]) == piecewise_constant(clocks=[0.0, 3.0, 5, 6], values=[3, 4, 5.0])
    @test piecewise_constant(clocks=[0, 3, 5, 6], values=[3, 4, 8]) != piecewise_constant(clocks=[0.0, 3.0, 5, 6], values=[3, 4, 5.0])
end

@testset "waveform equality" begin
    @test Waveform(piecewise_linear(clocks=[0, 3, 5], values=[3, 4, 5]), 10) == Waveform(piecewise_linear(clocks=[0, 3.0, 5], values=[3, 4, 5]), 10.0)
    @test Waveform(piecewise_linear(clocks=[0, 3, 5], values=[3, 4, 5]), 10) != Waveform(piecewise_constant(clocks=[0, 3, 5, 6], values=[3, 4, 8]), 10.0)
end


@testset "append" begin
    duration=1.0
    wf = constant(;duration=duration,value=1.0)
    max_slope = 10
    t_begin = 1.0/max_slope
    t_end = duration - 1.0/max_slope
    begin_value = 0.0
    end_value = 0.0
    target_wf = piecewise_linear(;clocks=[0.0,t_begin,t_end,1],values=[0.0,1.0,1.0,0.0])

    mid_wf = wf[t_begin..t_end]

    start_wf = linear_ramp(;
        duration=t_begin,
        start_value=begin_value,
        stop_value=wf(t_begin)
    )

    end_wf = linear_ramp(;
        duration=duration-t_end,
        start_value=wf(t_end),
        stop_value=end_value
    )
    new_wf = append(start_wf,mid_wf,end_wf)

    @test new_wf ≈ target_wf

    # check floating point error
    wf = append(
        piecewise_linear(clocks=[0.0,0.2],values=[0.0,0.0]),
        piecewise_linear(clocks=[0.0,0.16], values = [0.1,0.1]),
        piecewise_linear(clocks=[0.0,0.2], values=[0.0,0.0])
    );

    @test wf(0.56) ≈ 0.0
end
