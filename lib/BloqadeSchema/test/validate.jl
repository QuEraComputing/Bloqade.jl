using Test 
using Unitful
using BloqadeWaveforms:
    Waveform,
    piecewise_linear,
    piecewise_constant
using BloqadeSchema:
    TaskCapabilities,
    check_resolution,
    message,
    get_rydberg_capabilities,
    get_device_capabilities,
    set_resolution,
    check_durations,
    validate_Ω,
    validate_Δ,
    validate_ϕ,
    validate_lattice,
    validate
using BloqadeExpr:
    rydberg_h
using BloqadeLattices

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
    @test message(!=) == "is not equal to the"

end

@testset "validate_lattice" begin
    # example for each error
    dc = get_device_capabilities()


    # maximum qubits
    nqubits_max = dc.lattice.number_qubits_max
    atoms = generate_sites(SquareLattice(),17,16,scale=4)
    nqubits = length(atoms)
    violations = validate_lattice(atoms,dc)

    @test violations == Set(["$nqubits qubits $(message(>)) of $nqubits_max qubits"])

    # position resolution
    position_resolution = dc.lattice.geometry.position_resolution
    position = (9*position_resolution/10,0)
    atoms = [position]
    violations = validate_lattice(atoms,dc)
    @test violations == Set(["1th atom position $position not consistent with position resolution: $position_resolution"])

    # maximum width
    max_width = dc.lattice.area.width
    width = max_width+1
    atoms = [(0,0),(width,0)]
    violations = validate_lattice(atoms,dc)
    @test violations == Set(["total width $width μm $(message(>)) value of $max_width μm"])


    # maximum height
    max_height = dc.lattice.area.height
    height = max_height+1
    atoms = [(0,0),(0,height)]
    violations = validate_lattice(atoms,dc)
    @test violations == Set(["total height $height μm $(message(>)) value of $max_height μm"])

    # radial distance
    ## NOTE: minimum supported radial distance is equal to the minimum supported vertical distance
    radial_spacing_min = dc.lattice.geometry.spacing_radial_min
    vertical_resolution = dc.lattice.geometry.spacing_vertical_min
    radial_spacing = 0.9 * radial_spacing_min
    v_1 = (0,0)
    v_2 = (0,radial_spacing)
    atoms = [v_1,v_2]
    violations = validate_lattice(atoms,dc)
    @test violations == Set([
        "positions {1 => $v_1, 2 => $v_2} violate y minimum value of $vertical_resolution μs",
        "positions 1 => $v_1 and 2 => $v_2 are a distance of $radial_spacing μm apart which is $(message(<)) value of $radial_spacing_min μm",
    ])


    v_1 = (0,0)
    v_2 = (5,2*vertical_resolution-0.1)
    v_3 = (10,vertical_resolution-0.1)
    v_4 = (15,0)
    atoms = [v_1,v_2,v_3,v_4]
    error_sites = join([site=>atoms[site] for site in [1,3,4]],", ")
    violations = validate_lattice(atoms,dc)
        
    @test violations == Set(["positions {$(error_sites)} violate y minimum value of $vertical_resolution μs"])


end

@testset "validate_Ω" begin

    rydberg_capabilities = get_rydberg_capabilities()

    max_time = rydberg_capabilities.Ω.max_time
    # waveform duration cannot exceed maximum supported time
    Ω = piecewise_linear(;clocks = [0.0,1.0,2.0,3.0,max_time + 0.1], values = fill(0.0, 5))
    @test validate_Ω(Ω, rydberg_capabilities.Ω) == Set([
        "Ω(t) duration with value $(max_time + 0.1) μs exceeds maximum value of $max_time μs"
    ])

    # waveform cannot have duration smaller than minimum allowable time step
    # can give two warnings: waveform duration is shorter than supported minimum time step
    #                        AND minimum time step is shorter than supported minimum time step

    min_time_step = rydberg_capabilities.Ω.min_time_step
    Ω = piecewise_linear(clocks = [0.0, min_time_step - 0.01], values = [0.0, 0.0])
    @test validate_Ω(Ω, rydberg_capabilities.Ω)== Set([
        "Ω(t) duration with value $(min_time_step - 0.01) μs below minimum value of $min_time_step μs",
        "Ω(t) minimum step with value $(min_time_step - 0.01) μs below minimum value of $min_time_step μs"
    ])

    # waveform minimum time step cannot be smaller than supported minimum (0.01 microseconds)
    Ω = piecewise_linear(clocks = [0.0, 0.01, 0.011, 0.1], values=fill(0.0, 4))
    @test validate_Ω(Ω, rydberg_capabilities.Ω) == Set([
        "Ω(t) minimum step with value 0.001 μs below minimum value of $min_time_step μs"
    ])

    # initial value must be 0.0 radians
    Ω = piecewise_linear(;clocks=[0.0,0.1,1.0,2.0], values=[0.1,0.5,0.9,0.0])
    @test validate_Ω(Ω, rydberg_capabilities.Ω) == Set([
        "Ω(t) start value with value 0.1 rad⋅MHz is not equal to the value of 0.0 rad⋅MHz"
    ])

    # end value must be 0.0 radians
    Ω = piecewise_linear(;clocks=[0.0,0.1,1.0,2.0], values=[0.0,0.5,0.9,1.0])
    @test validate_Ω(Ω, rydberg_capabilities.Ω) == Set([
        "Ω(t) end value with value 1.0 rad⋅MHz is not equal to the value of 0.0 rad⋅MHz"
    ])

    supported_max_slope = rydberg_capabilities.Ω.max_slope
    end_value = set_resolution(supported_max_slope*1.01 * 0.05, rydberg_capabilities.Ω.value_resolution)
    slope = end_value/0.05
    Ω = piecewise_linear(clocks = [0.0, 0.05, 1.0, 4.0], values = [0.0, end_value, end_value, 0])
    @test validate_Ω(Ω, rydberg_capabilities.Ω) == Set([
        "Ω(t) maximum slope with value $(slope) rad⋅MHz/μs exceeds maximum value of $supported_max_slope rad⋅MHz/μs"
    ])

    supported_min_value = rydberg_capabilities.Ω.min_value
    Ω = piecewise_linear(clocks = [0.0, 2.0, 4.0], values = [0.0, supported_min_value-0.1,0.0])
    @test validate_Ω(Ω, rydberg_capabilities.Ω) == Set([
        "Ω(t) minimum value with value $(supported_min_value-0.1) rad⋅MHz below minimum value of $supported_min_value rad⋅MHz"
    ])

    supported_max_value = rydberg_capabilities.Ω.max_value
    # cannot go below hardware supported minimum value
    Ω = piecewise_linear(;clocks=[0.0,2.0,4.0], values=[0.0, supported_max_value+0.1, 0.0])
    @test validate_Ω(Ω, rydberg_capabilities.Ω) == Set([
        "Ω(t) maximum value with value $(supported_max_value+0.1) rad⋅MHz exceeds maximum value of $supported_max_value rad⋅MHz"
    ])

    # time values must be integer multiple of resolution (0.001 μs)
    Ω = piecewise_linear(;clocks=[0.0,1.0,2.0,3.0,3.9995], values=fill(0.0,5))
    @test validate_Ω(Ω, rydberg_capabilities.Ω) == Set([
        "Ω(t) clock 3.9995 μs is not consistent with resolution 0.001 μs."
    ])

    # waveform values must be integer multiple of resolution (5.0e-7 rad)
    resolution = rydberg_capabilities.Ω.value_resolution
    value = round(11.4*resolution;sigdigits=14)
    Ω = piecewise_linear(;clocks=[0.0, 0.1, 0.2, 0.3], values=[0.0, 10*resolution, value, 0.0])
    @test validate_Ω(Ω, rydberg_capabilities.Ω) == Set([
        "Ω(t) value $value rad is not consistent with resolution $resolution rad⋅MHz."
    ])

end

@testset "validate_Δ" begin

    rydberg_capabilities = get_rydberg_capabilities()

    # waveform duration cannot exceed maximum supported time
    Δ = piecewise_linear(;clocks = [0.0,1.0,2.0,3.0,4.1], values = fill(0.0, 5))
    @test validate_Δ(Δ, rydberg_capabilities.Δ) == Set([
        "Δ(t) duration with value 4.1 μs exceeds maximum value of 4.0 μs"
    ])

    # waveform cannot have duration smaller than minimum allowable time step
    # can give two warnings: waveform duration is shorter than supported minimum time step
    #                        AND minimum time step is shorter than supported minimum time step
    min_time_step = rydberg_capabilities.Δ.min_time_step
    Δ = piecewise_linear(clocks = [0.0,0.009], values = [0.0, 0.0])
    @test validate_Δ(Δ, rydberg_capabilities.Δ)== Set([
        "Δ(t) duration with value $(Δ.duration) μs below minimum value of $min_time_step μs",
        "Δ(t) minimum step with value $(Δ.duration) μs below minimum value of $min_time_step μs"
    ])

    # waveform minimum time step cannot be smaller than supported minimum (0.01 microseconds)
    Δ = piecewise_linear(clocks = [0.0, 0.01, 0.011, 0.1], values=fill(0.0, 4))
    @test validate_Δ(Δ, rydberg_capabilities.Δ) == Set([
        "Δ(t) minimum step with value 0.001 μs below minimum value of $min_time_step μs"
    ])
    
    # violating slope constraint also means violating maximum value constraint
    supported_max_slope = rydberg_capabilities.Δ.max_slope
    end_value = set_resolution(supported_max_slope*1.01 * 0.05, rydberg_capabilities.Δ.value_resolution)
    slope = end_value/0.05
    Δ = piecewise_linear(clocks = [0.0, 0.05, 1.0, 4.0], values = [0.0, end_value, end_value, 0.0])
    @test validate_Δ(Δ, rydberg_capabilities.Δ) == Set([
        "Δ(t) maximum value with value $end_value rad⋅MHz exceeds maximum value of $(rydberg_capabilities.Δ.max_value) rad⋅MHz",
        "Δ(t) maximum slope with value $slope rad⋅MHz/μs exceeds maximum value of $supported_max_slope rad⋅MHz/μs"
    ])

    supported_min_value = rydberg_capabilities.Δ.min_value
    Δ = piecewise_linear(clocks = [0.0, 4.0], values = [0.0, supported_min_value-0.1])
    @test validate_Δ(Δ, rydberg_capabilities.Δ) == Set([
        "Δ(t) minimum value with value $(supported_min_value-0.1) rad⋅MHz below minimum value of $supported_min_value rad⋅MHz"
    ])

    supported_max_value = rydberg_capabilities.Δ.max_value
    # cannot go beyond hardware supported maximum Δ value
    Δ = piecewise_linear(;clocks=[0.0,4.0], values=[0.0, supported_max_value+0.1])
    @test validate_Δ(Δ, rydberg_capabilities.Δ) == Set([
        "Δ(t) maximum value with value $(supported_max_value+0.1) rad⋅MHz exceeds maximum value of $supported_max_value rad⋅MHz"
    ])

    # time values must be integer multiple of resolution (0.001 μs)
    Δ = piecewise_linear(;clocks=[0.0,1.0,2.0,3.0,3.9995], values=fill(0.0,5))
    @test validate_Δ(Δ, rydberg_capabilities.Δ) == Set([
        "Δ(t) clock 3.9995 μs is not consistent with resolution 0.001 μs."
    ])

    # waveform values must be integer multiple of resolution (5.0e-7 rad)
    resolution = rydberg_capabilities.Δ.value_resolution
    value = round(11.4*resolution;sigdigits=14)
    Δ = piecewise_linear(;clocks=[0.0, 0.1, 0.2, 0.3], values=[0.0, 10*resolution, value, 0.0])
    @test validate_Δ(Δ, rydberg_capabilities.Δ) == Set([
        "Δ(t) value $value rad is not consistent with resolution $resolution rad⋅MHz."
    ])
end

@testset "validate_ϕ" begin

    rydberg_capabilities = get_rydberg_capabilities()

    # waveform duration cannot exceed maximum supported time
    ϕ = piecewise_constant(;clocks = [0.0,1.0,2.0,3.0,4.1], values = fill(0.0, 4))
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) duration with value 4.1 μs exceeds maximum value of 4.0 μs"
    ])

    # waveform cannot have duration smaller than minimum allowable time step
    # can give two warnings: waveform duration is shorter than supported minimum time step
    #                        AND minimum time step is shorter than supported minimum time step
    min_time_step = rydberg_capabilities.ϕ.min_time_step
    ϕ = piecewise_constant(clocks = [0.0, min_time_step - 0.01], values = [0.0])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ)== Set([
        "ϕ(t) duration with value $(ϕ.duration) μs below minimum value of $min_time_step μs",
        "ϕ(t) minimum step with value $(ϕ.duration) μs below minimum value of $min_time_step μs"
    ])

    # waveform minimum time step cannot be smaller than supported minimum (0.01 microseconds)
    ϕ = piecewise_constant(clocks = [0.0, 0.01, 0.011, 0.1], values=fill(0.0, 3))
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) minimum step with value 0.001 μs below minimum value of $min_time_step μs"
    ])

    supported_min_value = rydberg_capabilities.ϕ.min_value
    ϕ = piecewise_constant(clocks = [0.0, 2.0, 4.0], values = [0.0, supported_min_value-0.1])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) minimum value with value $(supported_min_value-0.1) rad below minimum value of $supported_min_value rad"
    ])

    supported_max_value = rydberg_capabilities.ϕ.max_value
    # cannot go below hardware supported minimum value
    ϕ = piecewise_constant(;clocks=[0.0, 2.0,4.0], values=[0.0, supported_max_value+0.1])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) maximum value with value $(supported_max_value+0.1) rad exceeds maximum value of $supported_max_value rad"
    ])

    # initial value must be 0.0 radians
    ϕ = piecewise_constant(;clocks=[0.0,0.1,1.0,2.0], values=[0.1,0.5,0.9])
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) start value with value 0.1 rad is not equal to the value of 0.0 rad"
    ])

    # time values must be integer multiple of resolution (0.001 μs)
    ϕ = piecewise_constant(;clocks=[0.0,1.0,2.0,3.0,3.9995], values=fill(0.0,4))
    @test validate_ϕ(ϕ, rydberg_capabilities.ϕ) == Set([
        "ϕ(t) clock 3.9995 μs is not consistent with resolution 0.001 μs."
    ])

    # waveform values must be integer multiple of resolution (5.0e-7 rad)
    resolution = rydberg_capabilities.ϕ.value_resolution
    value = 11.4*resolution
    ϕ = piecewise_constant(;clocks=[0.0, 0.1, 0.2, 0.3], values=[0.0, 10*resolution, value])
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
    ϕ = piecewise_constant(;clocks=Float64[0,1,2,3],values=Float64[0,1,-1])
    atoms = 5.0 * [i for i in 1:10]
    H = rydberg_h(atoms, Ω=Ω, Δ=Δ, ϕ=ϕ)
    violations = validate(H)
    @test isempty(violations)
    
    show(stdout, MIME"text/plain"(),violations)

end