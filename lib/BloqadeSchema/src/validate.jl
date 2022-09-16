

error_or_warn(warn::Bool,msg::String) = (warn ? @warn(msg) : error(msg))

# checks of x is integer multiple of res
# if it is return false otherwise return true
function check_resolution(res::Float64,x::Float64)
    x == 0 && return false

    x_div_res = (abs(x) / res)
    return !isapprox(round(x_div_res),x_div_res)
end

# TODO: implement this test based on validation in TaskManager
function validate_lattice(positions,warn::Bool,::DeviceCapabilities) end

# Waveform Validations all waveforms must be piecewise linear

# helper functions to generate messages
message(::typeof(>)) = "exceeds maximum"
message(::typeof(<)) = "below minimum"
message(::typeof(!=)) = "is not equal to"

# TODO: add explicit typing, wf::Waveform{BloqadeWaveforms.PiecewiseLinear{T,I},T} where {T<:Real,I}
function validate_Ω(wf,warn,expected)
    
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
        ("minimum step",min_time_step,>,max_time,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad⋅MHz/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad⋅MHz"),
        ("maximum value",max_value,>,expected.max_value,"rad⋅MHz"),
        ("start value",start_value,!=,0.0,"rad⋅MHz"),
        ("end value",end_value,!=,0.0,"rad⋅MHz"),
    ]
    #("duration", max_time, < expected.min_time_step, "μs")

    for (name,given,op,expected,units) in tests
        op(given,expected) && error_or_warn(warn,"Ω(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    map(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && error_or_warn(warn,"Ω(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    map(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && error_or_warn(warn,"Ω(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad⋅MHz.")
    end

end

function validate_Δ(wf,warn,expected)
    
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

    for (name,given,op,expected,units) in tests
        op(given,expected) && error_or_warn(warn,"Δ(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    map(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && error_or_warn(warn,"Δ(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    map(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && error_or_warn(warn,"Δ(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad⋅MHz.")
    end

end

function validate_ϕ(wf,warn,expected)

    max_time = wf.duration
    max_value = maximum(wf.f.values)
    min_value = minimum(wf.f.values)
    min_time_step  = round(minimum(diff(wf.f.clocks));sigdigits=14)
    max_slope = round(maximum(abs.(diff(wf.f.values))./diff(wf.f.clocks));sigdigits=14)
    start_value = wf.f.values[1]


    tests = [
        ("duration",max_time,>,expected.max_time,"μs"),
        ("duration",max_time,<,expected.min_time_step,"μs"),
        ("minimum step",min_time_step,<,expected.min_time_step,"μs"),
        ("minimum step",min_time_step,>,max_time,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad"),
        ("maximum value",max_value,>,expected.max_value,"rad"),
        ("start value",start_value,!=,0.0,"rad"),
    ]

    for (name,given,op,expected,units) in tests
        op(given,expected) && error_or_warn(warn,"ϕ(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    map(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && error_or_warn(warn,"ϕ(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    map(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && error_or_warn(warn,"ϕ(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad.")
    end
end

function validate_δ(wf,Δi,warn,expected)

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

    for (name,given,op,expected,units) in tests
        op(given,expected) && error_or_warn(warn,"δ(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    map(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && error_or_warn(warn,"δ(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    map(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && error_or_warn(warn,"δ(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad.")
    end

    map(Δi) do value
        check_resolution(expected.local_mask_resolution,value) && error_or_warn(warn,"Δi value $value  is not consistent with resolution $(expected.local_mask_resolution).")
    end

end

# public API exposed here
function validate(H::BloqadeExpr.RydbergHamiltonian;warn::Bool=false,device_capabilities::DeviceCapabilities=get_device_capabilities())
    (atoms,ϕ,Ω,Δ,δ,Δi) = schema_parse_rydberg_fields(H)

    validate_lattice(H.rydberg_term.atoms,warn,device_capabilities)
    rydberg_capabilities = get_rydberg_capabilities(;device_capabilities=device_capabilities)

    validate_ϕ(ϕ,warn,rydberg_capabilities.ϕ)
    validate_Ω(Ω,warn,rydberg_capabilities.Ω)
    validate_Δ(Δ,warn,rydberg_capabilities.Δ)

    if !isnothing(δ)
        validate_δ(δ,Δi,warn,rydberg_capabilities.δ)
    end

end
