using Test
using BloqadeWaveforms
using QuadGK: quadgk

slope_msg = "Requested tolerance for interpolation violates the slope constraint."
step_msg = "Requested tolerance for interpolation violates the step size constraint."
constraint_msg = "Interpolation requires either a tolerance constraint or a slope and step constraint."
warn_msg = "negative tolerance provided, taking absolute value."

function benchmark_func() 
    wf = Waveform(t->t^2,1)
    new_wf = piecewise_constant_interpolate(wf; atol=0,min_step=1e-3)
end

@testset "piecewise constant" begin

    @testset "constraint terminated" begin
        wf = Waveform(t->t^2,1)
        new_wf = piecewise_constant_interpolate(wf; atol=0,min_step=1e-3)
        @test norm(wf - new_wf) < 1e-1
    end

    @testset "general waveform" begin

        f_list = [(t->t^2,4.0),(t->sin(t),2π)]
        atol=1e-2
        for (function_number, (f,duration)) in enumerate(f_list)
            wf = Waveform(f,duration)
            new_wf = piecewise_constant_interpolate(wf;atol=atol)
    
            @test isapprox(wf,new_wf,atol=atol)
        end
    
        wf = Waveform(t->t^2,2)
    
        @test_logs (:warn,warn_msg) piecewise_linear_interpolate(wf;atol=-1e-5)
        @test_throws ErrorException piecewise_linear_interpolate(wf;min_step = 0.1)
        @test_throws ErrorException piecewise_linear_interpolate(wf;atol=0)
    
    end

    @testset "piecewise constant waveforms" begin
        wf = piecewise_constant(;clocks = [0.0, 2.0, 3.0, 4.0], values = [0.0, 2.0, 1.0])
        new_wf = piecewise_constant_interpolate(wf)
        @test wf == new_wf

        step_msg = "Waveform step smaller than constraint."
        @test_logs (:warn,warn_msg) piecewise_constant_interpolate(wf;atol=-1e-5)
        @test_throws ErrorException piecewise_constant_interpolate(wf;atol=0)
        @test_throws ErrorException piecewise_constant_interpolate(wf;min_step=3.0)
    end

end

@testset "piecewise linear waveforms" begin
    wf = piecewise_linear(clocks = [0.0, 2.0, 3.0, 4.0], values = [0.0, 3.0, 1.1, 2.2])
    new_wf = piecewise_linear_interpolate(wf)
    @test wf == new_wf
    
    slope_msg = "Waveform slope larger than constraint."
    step_msg = "Waveform step smaller than constraint."
    # test_log instead of test_warn for julia 1.6
    @test_logs (:warn,warn_msg) piecewise_linear_interpolate(wf;atol=-1e-5)
    @test_throws ErrorException piecewise_linear_interpolate(wf;atol=0)
    @test_throws ErrorException piecewise_linear_interpolate(wf;min_step = 10.0)
    @test_throws ErrorException piecewise_linear_interpolate(wf;max_slope = 0.1)

end

@testset "piecewise constant waveforms" begin
    wf = piecewise_constant(;clocks = [0.0, 2.0, 3.0, 4.0], values = [0.0, 2.0, 1.0])

    pwc_step_msg = "Distance between steps in piecewise constant waveform are too small to convert to piecewise linear."

    new_wf = piecewise_linear_interpolate(wf,atol=1e-3) # no constraints
    @test new_wf == piecewise_linear(
        clocks=[0.0,1.9995,2.0005,2.999,3.001,4.0],
        values=[0.0,0.0,2.0,2.0,1.0,1.0]
        )

    @test_logs (:warn,warn_msg) piecewise_linear_interpolate(wf;atol=-1e-5)
    @test_throws ErrorException piecewise_linear_interpolate(wf;atol=0)
    @test_throws ErrorException piecewise_linear_interpolate(wf;min_step = 1.0)
    @test_throws ErrorException piecewise_linear_interpolate(wf;max_slope = 1.0)


    wf = piecewise_constant(;clocks = [0.0, 1.0, 2.0], values = [0.0, 20.0])
    @test_throws ErrorException piecewise_linear_interpolate(wf;atol=0,max_slope=5.0,min_step=0.1)

end

@testset "general waveform" begin
    f_list = [(t->t^2,10.0),(t->sin(t),2π),(t->t*sin(t^2),10.0),(t->sqrt(t)*sign(sin(t)),2π)]
    atol=1e-3
    for (f,duration) in f_list
        wf = Waveform(f,duration)
        new_wf = piecewise_linear_interpolate(wf;atol=atol)


        @test isapprox(wf,new_wf,atol=atol)
    end

    wf = Waveform(t->t^2,2)

    @test_logs (:warn,warn_msg) piecewise_linear_interpolate(wf;atol=-1e-5)
    @test_throws ErrorException piecewise_linear_interpolate(wf;max_slope = 2.0)
    @test_throws ErrorException piecewise_linear_interpolate(wf;min_step = 0.1)
    @test_throws ErrorException piecewise_linear_interpolate(wf;atol=0)

end

@testset "constraint terminated" begin
    wf = Waveform(t->t^2,1)
    new_wf = piecewise_linear_interpolate(wf,atol=0,max_slope=100,min_step=1e-5)
    @test norm(wf - new_wf) < 1e-3
end