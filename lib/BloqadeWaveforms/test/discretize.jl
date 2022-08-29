using Test
using BloqadeWaveforms
using QuadGK: quadgk

@testset "piecewise linear waveforms" begin
    wf = piecewise_linear(clocks = [0.0, 2.0, 3.0, 4.0], values = [0.0, 3.0, 1.1, 2.2])
    new_wf = discretize(wf)
    @test wf == new_wf

    @test_throws ErrorException discretize(wf;max_value=1.0)
    @test_throws ErrorException discretize(wf;min_step = 10.0)
    @test_throws ErrorException discretize(wf;max_slope = 0.1)

end

@testset "piecewise constant waveforms" begin
    wf = piecewise_constant(clocks = [0.0, 2.0, 3.0, 4.0], values = [0.0, 2.0, 1.0])

    new_wf = discretize(wf,tol=1e-3) # no constraints
    @test new_wf == piecewise_linear(
        clocks=[0.0,1.9995,2.0005,2.999,3.001,4.0],
        values=[0.0,0.0,2.0,2.0,1.0,1.0]
        )

        @test_throws ErrorException discretize(wf;max_value=0.5)
        @test_throws ErrorException discretize(wf;min_step = 1.0)
        @test_throws ErrorException discretize(wf;max_slope = 1.0)

end

@testset "general waveform" begin
    f_list = [(t->t^2,10.0),(t->sin(t),2π),(t->t*sin(t^2),10.0),(t->sqrt(t)*sign(sin(t)),2π)]
    tol=1e-3
    for (f,duration) in f_list
        wf = Waveform(f,duration)
        new_wf = discretize(wf;tol=tol)


        @test isapprox(wf,new_wf,atol=tol,rtol=0)
    end

    wf = Waveform(t->t^2,2)

    @test_throws ErrorException discretize(wf;max_slope = 2.0)
    @test_throws ErrorException discretize(wf;min_step = 0.1)
    @test_throws ErrorException discretize(wf;max_value=3)



end