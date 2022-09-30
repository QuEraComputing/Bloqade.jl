


# checks of x is integer multiple of res
# if it is return false otherwise return true
function check_resolution(res::Float64,x::Float64)
    x == 0 && return false

    x_div_res = (abs(x) / res)
    return !isapprox(round(x_div_res),x_div_res)
end

# TODO: implement this test based on validation in TaskManager
validate_lattice(positions,::DeviceCapabilities) = Set([])

# Waveform Validations all waveforms must be piecewise linear

# helper functions to generate violations
message(::typeof(>)) = "exceeds maximum"
message(::typeof(<)) = "below minimum"
message(::typeof(!=)) = "is not equal to the"

# TODO: add explicit typing, wf::Waveform{BloqadeWaveforms.PiecewiseLinear{T,I},T} where {T<:Real,I}
function validate_Ω(wf,expected)
    
    max_time = wf.duration
    max_value = maximum(wf.f.values)
    min_value = minimum(wf.f.values)
    min_time_step  = round(minimum(diff(wf.f.clocks));sigdigits=14)
    max_slope = round(maximum(abs.(diff(wf.f.values))./diff(wf.f.clocks));sigdigits=14)
    start_value = wf.f.values[1]
    end_value = wf.f.values[end]

    tests = [
        ("duration",max_time,>,expected.max_time,"μs"),
        ("duration",max_time,<,expected.min_time_step,"μs"),
        ("minimum step",min_time_step,<,expected.min_time_step,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad⋅MHz/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad⋅MHz"),
        ("maximum value",max_value,>,expected.max_value,"rad⋅MHz"),
        ("start value",start_value,!=,0.0,"rad⋅MHz"),
        ("end value",end_value,!=,0.0,"rad⋅MHz"),
    ]
    #("duration", max_time, < expected.min_time_step, "μs")
    violations = Set([])

    for (name,given,op,expected,units) in tests
        op(given,expected) && push!(violations,"Ω(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    foreach(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && push!(violations,"Ω(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    foreach(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && push!(violations,"Ω(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad⋅MHz.")
    end

    return violations

end

function validate_Δ(wf,expected)
    
    max_time = wf.duration
    max_value = maximum(wf.f.values)
    min_value = minimum(wf.f.values)
    min_time_step  = round(minimum(diff(wf.f.clocks));sigdigits=14)
    max_slope = round(maximum(abs.(diff(wf.f.values))./diff(wf.f.clocks));sigdigits=14)


    tests = [
        ("duration",max_time,>,expected.max_time,"μs"),
        ("duration",max_time,<,expected.min_time_step,"μs"),
        ("minimum step",min_time_step,<,expected.min_time_step,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad⋅MHz/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad⋅MHz"),
        ("maximum value",max_value,>,expected.max_value,"rad⋅MHz"),
    ]

    violations = Set([])

    for (name,given,op,expected,units) in tests
        op(given,expected) && push!(violations,"Δ(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    foreach(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && push!(violations,"Δ(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    foreach(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && push!(violations,"Δ(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad⋅MHz.")
    end

    return violations

end

function validate_ϕ(wf,expected)

    max_time = wf.duration
    max_value = maximum(wf.f.values)
    min_value = minimum(wf.f.values)
    min_time_step  = round(minimum(diff(wf.f.clocks));sigdigits=14)
    start_value = wf.f.values[1]


    tests = [
        ("duration",max_time,>,expected.max_time,"μs"),
        ("duration",max_time,<,expected.min_time_step,"μs"),
        ("minimum step",min_time_step,<,expected.min_time_step,"μs"),
        ("minimum value",min_value,<,expected.min_value,"rad"),
        ("maximum value",max_value,>,expected.max_value,"rad"),
        ("start value",start_value,!=,0.0,"rad"),
    ]

    violations = Set([])

    for (name,given,op,expected,units) in tests
        op(given,expected) && push!(violations,"ϕ(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    foreach(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && push!(violations,"ϕ(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    foreach(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && push!(violations,"ϕ(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad.")
    end

    return violations
end

function validate_δ(wf,Δi,expected)

    max_time = wf.duration
    max_value = maximum(wf.f.values)
    min_value = minimum(wf.f.values)
    min_time_step  = round(minimum(diff(wf.f.clocks));sigdigits=14)
    max_slope = round(maximum(abs.(diff(wf.f.values))./diff(wf.f.clocks));sigdigits=14)


    tests = [
        ("duration",max_time,>,expected.max_time,"μs"),
        ("duration",max_time,<,expected.min_time_step,"μs"),
        ("minimum step",min_time_step,<,expected.min_time_step,"μs"),
        ("minimum step",min_time_step,>,max_time,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad⋅MHz/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad⋅MHz"),
        ("maximum value",max_value,>,expected.max_value,"rad⋅MHz"),
    ]

    violations = Set([])

    for (name,given,op,expected,units) in tests
        op(given,expected) && push!(violations,"δ(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    foreach(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && push!(violations,"δ(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    foreach(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && push!(violations,"δ(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad.")
    end

    foreach(Δi) do value
        check_resolution(expected.local_mask_resolution,value) && push!(violations,"Δi value $value  is not consistent with resolution $(expected.local_mask_resolution).")
    end

    return violations
end

function check_durations(ϕ,Ω,Δ,δ)
    durations = Dict(:Δ=>Δ.duration,:Ω=>Ω.duration,:ϕ=>ϕ.duration)
    
    invalid = Set([])

    if !isnothing(δ)
        durations[:δ] = δ.duration
    end

    for (f1,d1) in durations
        for (f2,d2) in durations
           d1!=d2 && push!(invalid,Set([f1=>d1,f2=>d2]))
        end
    end

    violations = Set([])

    for ((f1,d1),(f2,d2)) in invalid
        if d2>d1
            push!(violations,"$f2(t) duration of $d2 μs is greater than $f1(t) duration of $d1 μs")
        else
            push!(violations,"$f1(t) duration of $d1 μs is greater than $f2(t) duration of $d2 μs")
        end
    end

    return violations
end

function validate_analog_params(atoms,ϕ,Ω,Δ,δ,Δi,device_capabilities::DeviceCapabilities)
    
    lattice_violations = validate_lattice(atoms,device_capabilities)
    rydberg_capabilities = get_rydberg_capabilities(;device_capabilities=device_capabilities)

    ϕ_violations = validate_ϕ(ϕ,rydberg_capabilities.ϕ)
    Ω_violations = validate_Ω(Ω,rydberg_capabilities.Ω)
    Δ_violations = validate_Δ(Δ,rydberg_capabilities.Δ)
    
    δ_violations = if !isnothing(δ)
        validate_δ(δ,Δi,rydberg_capabilities.δ)
    else
        Set([])
    end

    misc_violations = check_durations(ϕ,Ω,Δ,δ)

    return ValidationViolations(;lattice_violations=lattice_violations,
        misc_violations=misc_violations,
        Δ_violations=Δ_violations,
        δ_violations=δ_violations,
        Ω_violations=Ω_violations,
        ϕ_violations=ϕ_violations
    )
end

# public API exposed here
function validate(H::BloqadeExpr.RydbergHamiltonian;device_capabilities::DeviceCapabilities=get_device_capabilities())
    (atoms,ϕ,Ω,Δ,δ,Δi) = schema_parse_rydberg_fields(H)

    return validate_analog_params(atoms,ϕ,Ω,Δ,δ,Δi,device_capabilities)

end
