

error_or_warn(warn::Bool,msg::String) = (warn ? @warn(msg) : error(msg))

# checks of x is integer multiple of res
# if it is return false otherwise return true
function check_resolution(res::Float64,x::Float64)
    x == 0 && return false

    x_div_res = (abs(x) / res)
    return !isapprox(round(x_div_res),x_div_res)
end

# TODO: implement this test based on validation in TaskManager
validate_lattice(positions,warn::Bool,::DeviceCapabilities) = String[]

# Waveform Validations all waveforms must be piecewise linear

# helper functions to generate messages
message(::typeof(>)) = "exceeds maximum"
message(::typeof(<)) = "below minimum"
message(::typeof(!=)) = "is not equal to"

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
        ("minimum step",min_time_step,>,max_time,"μs"),
        ("maximum slope",max_slope,>,expected.max_slope,"rad⋅MHz/μs"),
        ("minimum value",min_value,<,expected.min_value,"rad⋅MHz"),
        ("maximum value",max_value,>,expected.max_value,"rad⋅MHz"),
        ("start value",start_value,!=,0.0,"rad⋅MHz"),
        ("end value",end_value,!=,0.0,"rad⋅MHz"),
    ]
    #("duration", max_time, < expected.min_time_step, "μs")
    messages = String[]

    for (name,given,op,expected,units) in tests
        op(given,expected) && push!(messages,"Ω(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    foreach(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && push!(messages,"Ω(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    foreach(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && push!(messages,"Ω(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad⋅MHz.")
    end

    return messages

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

    messages = String[]

    for (name,given,op,expected,units) in tests
        op(given,expected) && push!(messages,"Δ(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    foreach(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && push!(messages,"Δ(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    foreach(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && push!(messages,"Δ(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad⋅MHz.")
    end

    return messages

end

function validate_ϕ(wf,expected)

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

    messages = String[]

    for (name,given,op,expected,units) in tests
        op(given,expected) && push!(messages,"ϕ(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    foreach(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && push!(messages,"ϕ(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    foreach(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && push!(messages,"ϕ(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad.")
    end

    return messages
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

    messages = String[]

    for (name,given,op,expected,units) in tests
        op(given,expected) && push!(messages,"δ(t) $name with value $given $units $(message(op)) value of $expected $units")
    end

    foreach(wf.f.clocks) do clock
        check_resolution(expected.time_resolution,clock) && push!(messages,"δ(t) clock $clock μs is not consistent with resolution $(expected.time_resolution) μs.")
    end

    foreach(wf.f.values) do value
        check_resolution(expected.value_resolution,value) && push!(messages,"δ(t) value $value rad is not consistent with resolution $(expected.value_resolution) rad.")
    end

    foreach(Δi) do value
        check_resolution(expected.local_mask_resolution,value) && push!(messages,"Δi value $value  is not consistent with resolution $(expected.local_mask_resolution).")
    end

    return messages
end

function check_durations(ϕ,Ω,Δ,δ,warn)
    durations = Dict(:Δ=>Δ.duration,:Ω=>Ω.duration,:ϕ=>ϕ.duration)

    if !isnothing(δ)
        durations[:δ] = δ.duration
    end

    for (f1,d1) in durations
        for (f2,d2) in durations
           d1!=d2 && push!(messages,"$f1(t) duration of $d1 μs is not equal to $f2(t) duration of $d2 μs")
        end
    end
end

function validate_analog_fields(atoms,ϕ,Ω,Δ,δ,Δi,warn::Bool,device_capabilities::DeviceCapabilities)
    
    lattice_messages = validate_lattice(atoms,warn,device_capabilities)
    rydberg_capabilities = get_rydberg_capabilities(;device_capabilities=device_capabilities)

    Δ_messages = validate_ϕ(ϕ,rydberg_capabilities.ϕ)
    Ω_messages = validate_Ω(Ω,rydberg_capabilities.Ω)
    ϕ_messages = validate_Δ(Δ,rydberg_capabilities.Δ)
    
    δ_messages = if !isnothing(δ)
        validate_δ(δ,Δi,rydberg_capabilities.δ)
    end

    misc_messages = check_durations(ϕ,Ω,Δ,δ,warn)

    return 
end

# public API exposed here
function validate(H::BloqadeExpr.RydbergHamiltonian;warn::Bool=false,device_capabilities::DeviceCapabilities=get_device_capabilities())
    (atoms,ϕ,Ω,Δ,δ,Δi) = schema_parse_rydberg_fields(H)

    validate_analog_fields(atoms,ϕ,Ω,Δ,δ,Δi,warn,device_capabilities)

end
