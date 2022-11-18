using GarishPrint
using JSON

const PiecewiseLinearWaveform = Waveform{BloqadeWaveforms.PiecewiseLinear{T,I},T} where {T<:Real,I}
const PiecewiseConstantWaveform = Waveform{BloqadeWaveforms.PiecewiseConstant{T},T} where {T<:Real}


@option struct SchemaTranslationParams <: QuEraSchema
    n_shots::Int = 1
    device_capabilities::DeviceCapabilities = get_device_capabilities()
end

@Base.kwdef struct HardwareTransformInfo <: QuEraSchema
    ϕ_error
    Ω_error
    Δ_error
    Δ_mask
    mse_atoms
end


@option struct ShotOutput <: QuEraSchema
    shot_status_code::Int32
    pre_sequence::Vector{Int32}
    post_sequence::Vector{Int32}
end

"""
    struct TaskOutput <: QuEraSchema

The result of executing a `TaskSpecification` on the machine.

Output of [`execute`](@ref) function.

# Fields
- `task_status_code::Int`: Task Status
- `shot_outputs::Vector{ShotOutput}`: Contains pre- and post- shot 
sequence in binary of if atoms are in Rydberg/Ground state.
"""
@option struct TaskOutput <: QuEraSchema
    task_status_code::Int
    shot_outputs::Vector{ShotOutput}
end


Base.@kwdef struct ValidationViolations <: QuEraSchema
    lattice_violations::Set
    misc_violations::Set
    Δ_violations::Set
    Ω_violations::Set
    ϕ_violations::Set
    δ_violations::Set
end

struct ValidationException <: Exception 
    violations::ValidationViolations
end

Base.isempty(t::ValidationViolations) = (
        isempty(t.lattice_violations) && 
        isempty(t.misc_violations)    && 
        isempty(t.Δ_violations)       && 
        isempty(t.Ω_violations)       && 
        isempty(t.ϕ_violations)       && 
        isempty(t.δ_violations)
    )



function Base.show(io::IO, ::MIME"text/plain", t::ValidationViolations)
    violations_list = [
        t.lattice_violations, 
        t.misc_violations,
        t.ϕ_violations, 
        t.Ω_violations,
        t.Δ_violations, 
        t.δ_violations
    ]
    println(io,"The following validation violations occured:\n")
    nviolation = one(Int)
    foreach(violations_list) do violations
        for msg in violations
            println(io,nviolation,". ",msg)
            nviolation+=1
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", t::TaskSpecification)
    GarishPrint.pprint_struct(t)
end

function Base.show(io::IO, ::MIME"text/plain", t::TaskOutput)
    GarishPrint.pprint_struct(t)
end

function Base.showerror(io::IO, e::ValidationException)
    show(io,MIME"text/plain"(),e.violations)
end