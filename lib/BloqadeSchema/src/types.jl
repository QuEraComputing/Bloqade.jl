using GarishPrint
using JSON

const PiecewiseLinearWaveform = Waveform{BloqadeWaveforms.PiecewiseLinear{T,I},T} where {T<:Real,I}
const PiecewiseConstantWaveform = Waveform{BloqadeWaveforms.PiecewiseConstant{T},T} where {T<:Real}

"""
    struct SchemaTranslationParams

Used to specify number of times a Hamiltonian should be executed and the capabilities of the machine the Hamiltonian
will be executed on as an argument to functions that convert Hamiltonians to other formats.

See also [`to_schema`](@ref), [`to_dict`](@ref), [`to_json`](@ref)

# Fields
- `nshots::Int = 1`:  The number of times a Hamiltonian should be executed
- `device_capabilities::DeviceCapabilities = get_device_capabilities()`: capabilities of the machine the Hamiltonian will be executed on
"""
@option struct SchemaTranslationParams <: QuEraSchema
    n_shots::Int = 1
    device_capabilities::DeviceCapabilities = get_device_capabilities()
end

"""
    struct HardwareTransformInfo <: QuEraSchema

Contains the calculated differences (error) betwen the original and transformed waveforms and atoms from
invoking [`hardware_transform`](@ref) on a [`BloqadeExpr.RydbergHamiltonian`](@ref).

# Fields
- `ϕ_error`: Error between the original laser phase waveform (``A``) and transformed one (``B``) waveforms calculated as ``\\Vert A - B\\Vert_1``. 
- `Ω_error`: Error between the original Rabi frequency waveform (``A``) and transformed one (``B``) waveforms calculated as ``\\Vert A - B\\Vert_1``.
- `Δ_error`: Error between the global detuning waveform (``A``) and transformed one (``B``) waveforms calculated as ``\\Vert A - B\\Vert_1``.
- `Δ_mask`: Decoupling of local detuning field inferred from the detuning value specified in Δ.
!!! note "Local Detuning Support"
    Local Detunings are currently not supported by Bloqade but will be in future releases.
- `mse_atoms`: Mean Squared Error between original atom positions and transformed ones.

"""
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

"""
    struct ValidationViolations <: QuEraSchema

Stores violations of hardware constraints from the user-supplied
[`BloqadeExpr.RydbergHamiltonian`](@ref) as strings in sets.
This is returned by [`validate`](@ref) and [`to_schema`](@ref).

# Fields
- `lattice_violations::Set`: violations of hardware-supported lattice geometry
- `misc_violations::Set`: violations that do not fall into other categories (e.g. number of shots)
- `Δ_violations::Set`: violations of detuning waveform
- `Ω_violations::Set`: violations of Rabi frequency waveform
- `ϕ_violations::Set`: violations of Phase waveform
- `δ_violations::Set`: violations of local detuning waveforms
"""
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