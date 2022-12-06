

# helper functions to generate violations

message(::typeof(>)) = "exceeds maximum"
message(::typeof(<)) = "below minimum"
message(::typeof(!=)) = "is not equal to the"

# checks of x is integer multiple of res
# if it is return false otherwise return true
function check_resolution(res::Real,x::Real) 
    x == 0 && return false

    x_div_res = (abs(x) / res)
    return !isapprox(round(x_div_res),x_div_res)
end

function group_by_key(itr,key)
    groups = Dict{Any,Set}()

    foreach(enumerate(itr)) do (i,value)
        key_value = value[key]
        if haskey(groups,key_value)
            push!(groups[key_value],i)
        else
            groups[key_value] = Set([i])
        end
    end

    return groups
end

# TODO: implement this test based on validation in TaskManager
function validate_lattice(positions,dc::DeviceCapabilities)

    positions = map(positions) do position
        if length(position) == 1
            (position[1],0.0)
        else
            position
        end
    end

    violations = Set([])

    # valudate number of sites (qubits)
    nqubits = length(positions)
    nqubits_max = dc.lattice.number_qubits_max
    nqubits > nqubits_max && push!(violations,
        "$nqubits qubits $(message(>)) of $nqubits_max qubits"
    )

    # validate resolution
    position_resolution = dc.lattice.geometry.position_resolution
    foreach(enumerate(positions)) do (i,position)
        any(check_resolution.(position_resolution,position)) && push!(violations,
            "$(i)th atom position $position not consistent with position resolution: $position_resolution"
        )
    end

    # validate area

    xmin = minimum(pos->pos[1],positions)
    xmax = maximum(pos->pos[1],positions)
    ymin = minimum(pos->pos[2],positions)
    ymax = maximum(pos->pos[2],positions)

    width = xmax-xmin

    height = ymax-ymin
    max_width = dc.lattice.area.width
    max_height = dc.lattice.area.height
    width > max_width && push!(violations,
        "total width $width μm $(message(>)) value of $max_width μm"
    )
    height > max_height && push!(violations,
        "total height $height μm $(message(>)) value of $max_height μm"
    )

    # validate radial spacing
    radial_spacing_min = dc.lattice.geometry.spacing_radial_min
    for i in 1:nqubits, j in i+1:nqubits
        x_i,y_i = positions[i]
        x_j,y_j = positions[j]
        radial_spacing = hypot(x_i-x_j,y_i-y_j)
        atom_1 = i => (x_i,y_i)
        atom_2 = j => (x_j,y_j)
        radial_spacing < radial_spacing_min && push!(violations,
            "positions $atom_1 and $atom_2 are a distance of $radial_spacing μm apart which is $(message(<)) value of $radial_spacing_min μm"
        )
        
    end

    # validate vertical spacing

    y_groups = group_by_key(positions,2)
    y_values = collect(keys(y_groups))
    sort!(y_values)
    zip_with_next = zip(y_values[1:end-1],y_values[2:end])
    vertical_resolution = dc.lattice.geometry.spacing_vertical_min

    foreach(zip_with_next) do (y0,y1)
        if !(abs(y0-y1) > vertical_resolution || abs(y0-y1) ≈ vertical_resolution)
            sites = collect(union(y_groups[y0],y_groups[y1]))
            sorted_sites = sort(sites)
            error_sites = join([site=>positions[site] for site in sorted_sites],", ")
            push!(violations,
                "positions {$(error_sites)} violate y minimum value of $vertical_resolution μs"
            )

        end
    end

    return violations

end

# Waveform Validations all waveforms must be piecewise linear


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
"""
    validate(H::BloqadeExpr.RydbergHamiltonian;device_capabilities::DeviceCapabilities=get_device_capabilities())

Checks if `H` is capable of being represented as part of the internal schema as well as if it falls
in the capabilities of what the machine can do via `device_capabilities`.

Returns [`ValidationViolations`](@ref) with each field containing a set of strings indicating which
constraints were violated for which part of `H`.

Violations include:

## Waveform Type
* ϕ is not of type `PiecewiseConstantWaveform`
* Ω and Δ are not of type `PiecewiseLinearWaveform`
## Atom Position
* Number of qubits requested exceeds what is supported by the device
* Atom positions exceed position resolution supported by the device
* The total width/height of the atom arrangement exceeds what is supported by the device
* The radial spacing between atoms is smaller than what is supported by the device
* The vertical row spacing between atoms is smaller than what is supported by the device
## General Waveform Constraints (apply to Ω, Δ, ϕ)
* duration exceeds device supported duration
* duration is smaller than device supported minimum time step
* smallest time step is smaller than supported smallest time step
* value is smaller than smallest supported value
* value is larger than largest supported value
## Ω Waveform specific constraints
* Slope exceeds largest supported slope
* Start and end values are not equal to 0.0 rad⋅MHz
## Δ Waveform specific constraints
* Slope exceeds largest supported slope
## ϕ Waveform specific constraints
* start value is not equal to 0.0 rad⋅MHz

# Logs/Warnings/Exceptions

The following exceptions can be thrown:
* ϕ is not of type `PiecewiseConstantWaveform`
* Ω and Δ are not of type `PiecewiseLinearWaveform`

# Examples
```jldoctest; setup:=(using BloqadeExpr, BloqadeWaveforms, BloqadeLattices)
julia> Δ = Ω = ϕ = sinusoidal(duration=2, amplitude=1.3*π);

julia> h = rydberg_h(atom_positions; Ω=Ω,Δ=Δ,ϕ=ϕ)

julia> transformed_h, _ = transform(h); # transform returns error info

julia> validate(transformed_h) # constrained by default value of `device_capabilities` argument
The following validation violations occurred:

1. positions 2 => (2.0, 0.0) and 3 => (3.0, 0.0) are a distance of 1.0 μm apart which is below minimum value of 4.0 μm
2. positions 1 => (1.1, 0.0) and 2 => (2.0, 0.0) are a distance of 0.8999999999999999 μm apart which is below minimum value of 4.0 μm
3. positions 1 => (1.1, 0.0) and 3 => (3.0, 0.0) are a distance of 1.9 μm apart which is below minimum value of 4.0 μm
``` 
"""
function validate(H::BloqadeExpr.RydbergHamiltonian;device_capabilities::DeviceCapabilities=get_device_capabilities())
    (atoms,ϕ,Ω,Δ,δ,Δi) = schema_parse_rydberg_fields(H)

    return validate_analog_params(atoms,ϕ,Ω,Δ,δ,Δi,device_capabilities)

end
