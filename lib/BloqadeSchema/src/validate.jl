
error_or_warn(warn::Bool,msg::String) = (warn ? @warn(msg) : error(msg))

# checks of x is integer multiple of res
# if it is return false otherwise return true
function check_resolution(res::Float64,x::Float64)
    x == 0 && return true

    x_div_res = (abs(x) / res)
    return isapprox(round(x_div_res,sigdigits=15),x_div_res)
end

# TODO: implement this test based on validation in TaskManager
function validate_lattice(positions,warn::Bool,C::DeviceCapabilities) end

# Waveform Validations all waveforms must be piecewise linear

# helper functions to generate messages
message(::typeof(>)) = "exceeds maximum"
message(::typeof(<)) = "below minimum"
message(::typeof(==)) = "is not equal to"

function validate_Ω(wf,warn,expected)
    
    max_time = wf.duration
    max_value = maximum(wf.f.values)
    min_value = minimum(wf.f.values)
    min_step  = minimum(diff(wf.f.clocks))
    max_slope = maximum(abs.(diff(wf.f.values))./diff(wf.f.clocks))
    start_value = wf.f.values[1]
    end_value = wf.f.values[end]

    tests = [
        ("duration",max_time,>,expected.max_time,"μs"),
        ("minimum step",min_step,<,expected.min_step,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad⋅MHs/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad⋅MHs"),
        ("maximum value",max_value,>,expected.max_value,"rad⋅MHs"),
        ("start value",start_value,==,0.0,"rad⋅MHs"),
        ("end value",end_value,==,0.0,"rad⋅MHs"),
    ]

    for (name,given,op,expected,units) in tests
        op(given,expected) && error_or_warn(warn,"Ω(t) $name with value $given $units $message(op) value of $expected $units")
    end

    map(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && error_or_warn(warn,"Ω(t) clock $clock μs is not consistent with resolution $expected.time_resolution μs.")
    end

    map(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && error_or_warn(warn,"Ω(t) value $value rad⋅MHz is not consistent with resolution $expected.value_resolution rad⋅MHz.")
    end

end

function validate_Δ(wf,warn,expected)
    
    max_time = wf.duration
    max_value = maximum(wf.f.values)
    min_value = minimum(wf.f.values)
    min_step  = minimum(diff(wf.f.clocks))
    max_slope = maximum(abs.(diff(wf.f.values))./diff(Ω.f.clocks))


    tests = [
        ("duration",max_time,>,expected.max_time,"μs"),
        ("minimum step",min_step,<,expected.min_step,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad⋅MHz/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad⋅MHz"),
        ("maximum value",max_value,>,expected.max_value,"rad⋅MHz"),
    ]

    for (name,given,op,expected,units) in tests
        op(given,expected) && error_or_warn(warn,"Δ(t) $name with value $given $units $message(op) value of $expected $units")
    end

    map(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && error_or_warn(warn,"Δ(t) clock $clock μs is not consistent with resolution $expected.time_resolution μs.")
    end

    map(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && error_or_warn(warn,"Δ(t) value $value rad⋅MHz is not consistent with resolution $expected.value_resolution rad⋅MHz.")
    end

end

function validate_ϕ(wf,warn,expected)

    max_time = wf.duration
    max_value = maximum(wf.f.values)
    min_value = minimum(wf.f.values)
    min_step  = minimum(diff(wf.f.clocks))
    max_slope = maximum(abs.(diff(wf.f.values))./diff(Ω.f.clocks))


    tests = [
        ("duration",max_time,>,expected.max_time,"μs"),
        ("minimum step",min_step,<,expected.min_step,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad"),
        ("maximum value",max_value,>,expected.max_value,"rad"),
    ]

    for (name,given,op,expected,units) in tests
        op(given,expected) && error_or_warn(warn,"ϕ(t) $name with value $given $units $message(op) value of $expected $units")
    end

    map(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && error_or_warn(warn,"ϕ(t) clock $clock μs is not consistent with resolution $expected.time_resolution μs.")
    end

    map(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && error_or_warn(warn,"ϕ(t) value $value rad is not consistent with resolution $expected.value_resolution rad.")
    end
end

function validate_δ(wf,Δi,warn,expected)

    max_time = wf.duration
    max_value = maximum(wf.f.values)
    min_value = minimum(wf.f.values)
    min_step  = minimum(diff(wf.f.clocks))
    max_slope = maximum(abs.(diff(wf.f.values))./diff(Ω.f.clocks))


    tests = [
        ("duration",max_time,>,expected.max_time,"μs"),
        ("minimum step",min_step,<,expected.min_step,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad⋅MHz/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad⋅MHz"),
        ("maximum value",max_value,>,expected.max_value,"rad⋅MHz"),
    ]

    for (name,given,op,expected,units) in tests
        op(given,expected) && error_or_warn(warn,"δ(t) $name with value $given $units $message(op) value of $expected $units")
    end

    map(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && error_or_warn(warn,"δ(t) clock $clock μs is not consistent with resolution $expected.time_resolution μs.")
    end

    map(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && error_or_warn(warn,"δ(t) value $value rad⋅MHz is not consistent with resolution $expected.value_resolution rad.")
    end

    map(Δi) do value
        check_resolution(expected.local_mask_resolution,value) && error_or_warn(warn,"Δi value $value  is not consistent with resolution $expected.local_mask_resolution.")
    end

end

# public API exposed here
function validate(H::BloqadeExpr.RydbergHamiltonian;warn::Bool=false,device_capabilities::DeviceCapabilities=get_device_capabilities())
    (atoms,ϕ,Ω,Δ,δ,Δi) = schema_parse_rydberg_fields(H)

    validate_lattice(H.rydberg_term.atoms,warn,device_capabilities.lattice)
    rydberg_capabilities = get_rydberg_capabilities(;device_capabilities=device_capabilities)

    validate_ϕ(ϕ,warn,rydberg_capabilities.ϕ)
    validate_Ω(Ω,warn,rydberg_capabilities.Ω)
    validate_Δ(Δ,warn,rydberg_capabilities.Δ)
    validate_δ(δ,Δi,warn,rydberg_capabilities.δ)

end
