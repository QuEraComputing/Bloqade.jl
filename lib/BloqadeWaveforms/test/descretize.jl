using Test
using BloqadeWaveforms
using QuadGK: quadgk

@testset "piecewise linear waveforms" begin
    wf = piecewise_linear(clocks = [0.0, 2.0, 3.0, 4.0], values = [0.0, 3.0, 1.1, 2.2])
    new_wf = descretize_waveform(wf)
    @test wf == new_wf
     
    exception_max_value_test = false
    try
        descretize_waveform(wf;max_value=1.0)
    catch e
        if isa(e,ErrorException) && e.msg == "Waveform exceeds maximum value."
            exception_max_value_test = true
        end
    end

    exception_min_step_test = false
    try
        descretize_waveform(wf;min_step = 10.0)
    catch e
        if isa(e,ErrorException) && e.msg == "Waveform step smaller than constraint."
            exception_min_step_test = true
        end
    end

    exception_max_slope_test = false
    try
        descretize_waveform(wf;max_slope = 0.1)
    catch e
        if isa(e,ErrorException) && e.msg == "Waveform slope larger than constraint."
            exception_max_slope_test = true
        end
    end

    @test exception_max_value_test
    @test exception_min_step_test
    @test exception_max_slope_test

end

@testset "piecewise constant waveforms" begin
    wf = piecewise_constant(clocks = [0.0, 2.0, 3.0, 4.0], values = [0.0, 2.0, 1.0])

    new_wf = descretize_waveform(wf,tol=1e-3) # no constraints
    @test new_wf == piecewise_linear(
        clocks=[0.0,1.9995,2.0005,2.999,3.001,4.0],
        values=[0.0,0.0,2.0,2.0,1.0,1.0]
        )

        exception_max_value_test = false
        try
            descretize_waveform(wf;max_value=0.5)
        catch e
            if isa(e,ErrorException) && e.msg == "Waveform exceeds maximum value."
                exception_max_value_test = true
            end
        end
    
        exception_min_step_test = false
        try
            descretize_waveform(wf;min_step = 1.0)
        catch e
            if isa(e,ErrorException) && e.msg == "Descretization cannot obtain requested tolerance given the step size and slope constraints."
                exception_min_step_test = true
            end
        end

        exception_max_slope_test = false
        try
            descretize_waveform(wf;max_slope = 1.0)
        catch e
            if isa(e,ErrorException) && e.msg == "Descretization cannot obtain requested tolerance given the step size and slope constraints."
                exception_max_slope_test = true
            end
        end


        @test exception_max_value_test
        @test exception_min_step_test
        @test exception_max_slope_test

end

@testset "general waveform" begin
    f_list = [(t->t^2,10),(t->sin(t),2π),(t->t*sin(t^2),10),(t->sqrt(t)*sign(sin(t)),2π)]
    tol=1e-3
    for (f,duration) in f_list
        wf = Waveform(f,duration)
        new_wf = descretize_waveform(wf;tol=tol)

        error,_ = quadgk(t->abs.(wf(t).-new_wf(t)),0,wf.duration)

        @test error < tol
    end

    wf = Waveform(t->t^2,2)

    exception_max_value_test = false
    try
        descretize_waveform(wf;max_value=3)
    catch e
        if isa(e,ErrorException) && e.msg == "Waveform exceeds maximum value."
            exception_max_value_test = true
        end
    end

    exception_min_step_test = false
    try
        descretize_waveform(wf;min_step = 0.1)
    catch e
        if isa(e,ErrorException) && e.msg == "Descretization cannot obtain requested tolerance given the step size and slope constraints."
            exception_min_step_test = true
        end
    end

    exception_max_slope_test = false
    try
        descretize_waveform(wf;max_slope = 2.0)
    catch e
        if isa(e,ErrorException) && e.msg == "Descretization cannot obtain requested tolerance given the step size and slope constraints."
            exception_max_slope_test = true
        end
    end


    @test exception_max_value_test
    @test exception_min_step_test
    @test exception_max_slope_test

end