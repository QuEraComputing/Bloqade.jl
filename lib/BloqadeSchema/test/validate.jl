using Test 
using Unitful
using BloqadeWaveforms:
    Waveform,
    piecewise_linear
using BloqadeSchema:
    error_or_warn,
    check_resolution,
    message,
    get_rydberg_capabilities,
    get_device_capabilities,
    check_durations,
    validate_Ω,
    validate_Δ,
    validate_ϕ,
    validate
using BloqadeExpr:
    rydberg_h


@testset "check_resolution" begin
    # check if x is integer multiple of res,
    # return false if x is, otherwise return true
    @test check_resolution(0.001, 0.0) == false
    @test check_resolution(0.001, 10.0) == false
    @test check_resolution(0.037, 9.0)

end

@testset "message" begin

    @test message(>) == "exceeds maximum"
    @test message(<) == "below minimum"
    @test message(!=) == "is not equal to"

end

@testset "validate_Ω" begin

    rydberg_capabilities = get_rydberg_capabilities()
    # warn = false
    warn = false

    # waveform duration cannot exceed maximum supported time
    Ω = piecewise_linear(;clocks = [0.0,1.0,2.0,3.0,4.1], values = fill(0.0, 5))
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω) 

    # waveform cannot have duration smaller than minimum allowable time step
    # can give two warnings: waveform duration is shorter than supported minimum time step
    #                        AND minimum time step is shorter than supported minimum time step
    Ω = piecewise_linear(clocks = [0.0,0.009], values = [0.0, 0.0])
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω)

    # waveform minimum time step cannot be smaller than supported minimum (0.01 microseconds)
    Ω = piecewise_linear(clocks = [0.0, 0.009, 0.010, 0.1], values=fill(0.0, 4))
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω)

    # cannot have change in values larger than what system supports
    # 250.0 rad*MHz/μs
    Ω = piecewise_linear(;clocks=[0.0, 0.01, 0.02, 0.03], values=[0.0,25.0,25.0,0.0])
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω)

    # cannot go below hardware supported minimum value 
    supported_min_value = rydberg_capabilities.Ω.min_value
    Ω = piecewise_linear(;clocks=[0.0,0.1,1.0,2.0], values=[0.0,supported_min_value-0.01,1.0,0.0])
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω)

    # cannot exceed hardware supported maximum value
    supported_max_value = rydberg_capabilities.Ω.max_value
    Ω = piecewise_linear(;clocks=[0.0,0.1,1.0,2.0], values=[0.0,0.5,supported_max_value+0.01,0.0])
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω)

    # initial value must be 0.0
    Ω = piecewise_linear(;clocks=[0.0,0.1,1.0,2.0], values=[0.1,0.5,0.9,0.0])
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω)

    # final value must be 0.0 
    Ω = piecewise_linear(;clocks=[0.0,0.1,1.0,2.0], values=[0.0,0.5,0.9,0.1])
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω)

    # time values must be integer multiple of resolution (0.001 μs)
    Ω = piecewise_linear(;clocks=[0.0,1.0,2.0,3.0,3.9995], values=fill(0.0,5))
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω)

    # waveform values must be integer multiple of resolution
    # must be int multiples of 0.0004 rad*MHz
    Ω = piecewise_linear(;clocks=[0.0,0.1,0.2,0.3,0.4], values=[0.0, 0.0004, 0.0008, 0.0009, 0.0])
    @test validate_Ω(Ω, warn, rydberg_capabilities.Ω)
end

@testset "validate_Δ" begin

    rydberg_capabilities = get_rydberg_capabilities()
    warn = false

    # waveform duration cannot exceed maximum supported time
    Δ = piecewise_linear(;clocks = [0.0,1.0,2.0,3.0,4.1], values = fill(0.0, 5))
    @test validate_Δ(Δ, warn, rydberg_capabilities.Δ) 

    # waveform cannot have duration smaller than minimum allowable time step
    # can give two warnings: waveform duration is shorter than supported minimum time step
    #                        AND minimum time step is shorter than supported minimum time step
    Δ = piecewise_linear(clocks = [0.0,0.009], values = [0.0, 0.0])
    @test validate_Δ(Δ, warn, rydberg_capabilities.Δ)

    # waveform minimum time step cannot be smaller than supported minimum (0.01 microseconds)
    Δ = piecewise_linear(clocks = [0.0, 0.009, 0.010, 0.1], values=fill(0.0, 4))
    @test validate_Δ(Δ, warn, rydberg_capabilities.Δ)

    # cannot have change in values larger than what system supports (2500.0 rad*MHz/μs)
    # min_time_step of 0.01 μs
    supported_min_value = rydberg_capabilities.Δ.min_value
    supported_max_value = rydberg_capabilities.Δ.max_value
    Δ = piecewise_linear(clocks = [0.0, 0.01, 0.02, 0.03], values = [0.0, supported_min_value, supported_max_value, 0.0 ])
    @test validate_Δ(Δ, warn, rydberg_capabilities.Δ)

    # cannot go below hardware supported minimum value
    Δ = piecewise_linear(;clocks=[0.0,0.1,1.0,2.0], values=[supported_min_value-0.01,0.5,1.0,0.1])
    @test validate_Δ(Δ, warn, rydberg_capabilities.Δ)

    # cannot exceed hardware supported maximum value
    Δ = piecewise_linear(;clocks=[0.0,0.1,1.0,2.0], values=[0.0,0.5,supported_max_value+0.01,0.1])
    @test validate_Δ(Δ, warn, rydberg_capabilities.Δ)

    # time values must be integer multiple of resolution (0.001 μs)
    Δ = piecewise_linear(;clocks=[0.0,1.0,2.0,3.0,3.9995], values=fill(0.0,5))
    @test validate_Δ(Δ, warn, rydberg_capabilities.Δ)

    # waveform values must be integer multiple of resolution (2.0e-7 rad*MHz)
    resolution = rydberg_capabilities.Δ.value_resolution
    Δ = piecewise_linear(;clocks=[0.0, 0.1, 0.2, 0.3], values=[0.0, 10*resolution, 11.4*resolution, 0.0])
    @test validate_Δ(Δ, warn, rydberg_capabilities.Δ)

end

@testset "validate_ϕ" begin

    rydberg_capabilities = get_rydberg_capabilities()
    warn = false

    # waveform duration cannot exceed maximum supported time
    ϕ = piecewise_linear(;clocks = [0.0,1.0,2.0,3.0,4.1], values = fill(0.0, 5))
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) duration with value 4.1 μs exceeds maximum value of 4.0 μs"
    ])

    # waveform cannot have duration smaller than minimum allowable time step
    # can give two warnings: waveform duration is shorter than supported minimum time step
    #                        AND minimum time step is shorter than supported minimum time step
    ϕ = piecewise_linear(clocks = [0.0,0.009], values = [0.0, 0.0])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ)== Set([
        "ϕ(t) duration with value 0.009 μs below minimum value of 0.01 μs",
        "ϕ(t) minimum step with value 0.009 μs below minimum value of 0.01 μs"
    ])

    # waveform minimum time step cannot be smaller than supported minimum (0.01 microseconds)
    ϕ = piecewise_linear(clocks = [0.0, 0.01, 0.011, 0.1], values=fill(0.0, 4))
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) minimum step with value 0.001 μs below minimum value of 0.01 μs"
    ])

    supported_max_slope = rydberg_capabilities.ϕ.max_slope
    slope = round(supported_max_slope*1.01;sigdigits=10)
    end_value = slope * 0.1
    ϕ = piecewise_linear(clocks = [0.0, 0.1], values = [0.0, end_value])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) maximum slope with value $(supported_max_slope*1.01) rad/μs exceeds maximum value of $supported_max_slope rad/μs"
    ])

    supported_min_value = rydberg_capabilities.ϕ.min_value
    ϕ = piecewise_linear(clocks = [0.0, 4.0], values = [0.0, supported_min_value-0.1])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) minimum value with value $(supported_min_value-0.1) rad below minimum value of $supported_min_value rad"
    ])

    supported_max_value = rydberg_capabilities.ϕ.max_value
    # cannot go below hardware supported minimum value
    ϕ = piecewise_linear(;clocks=[0.0,4.0], values=[0.0, supported_max_value+0.1])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) maximum value with value $(supported_max_value+0.1) rad exceeds maximum value of $supported_max_value rad"
    ])

    # initial value must be 0.0 radians
    ϕ = piecewise_linear(;clocks=[0.0,0.1,1.0,2.0], values=[0.1,0.5,0.9,0.0])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) start value with value 0.1 rad is not equal to the value of 0.0 rad"
    ])

    # time values must be integer multiple of resolution (0.001 μs)
    ϕ = piecewise_linear(;clocks=[0.0,1.0,2.0,3.0,3.9995], values=fill(0.0,5))
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) clock 3.9995 μs is not consistent with resolution 0.001 μs."
    ])

    # waveform values must be integer multiple of resolution (5.0e-7 rad)
    resolution = rydberg_capabilities.ϕ.value_resolution
    value = 11.4*resolution
    ϕ = piecewise_linear(;clocks=[0.0, 0.1, 0.2, 0.3], values=[0.0, 10*resolution, value, 0.0])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) value $value rad is not consistent with resolution $resolution rad."
    ])

end

@testset "check_durations" begin
    δ = Waveform(t->t,3)
    Δ = Waveform(t->t,3)
    ϕ = Waveform(t->t,3)
    Ω = Waveform(t->t,4)


    @test check_durations(ϕ,Ω,Δ,δ) == Set([
        "Ω(t) duration of 4 μs is greater than δ(t) duration of 3 μs",
        "Ω(t) duration of 4 μs is greater than Δ(t) duration of 3 μs",
        "Ω(t) duration of 4 μs is greater than ϕ(t) duration of 3 μs",
    ])
    @test check_durations(ϕ,Ω,Δ,nothing) == Set([
        "Ω(t) duration of 4 μs is greater than Δ(t) duration of 3 μs",
        "Ω(t) duration of 4 μs is greater than ϕ(t) duration of 3 μs",
    ])
    
end

@testset "validate" begin

    Ω = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,1,0])
    Δ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[1,1,-1,-1])
    ϕ = piecewise_linear(;clocks=Float64[0,1,2,3],values=Float64[0,1,-1,-1])
    atoms = 5.0 * [i for i in 1:10]
    H = rydberg_h(atoms, Ω=Ω, Δ=Δ, ϕ=ϕ)
    @test_nowarn validate(H)

end