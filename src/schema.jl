struct ValidationError <: Exception
    msg::String
end

const PulseParameter = Maybe{Union{Float64, Vector{Float64}}}

function check_pulse_parameters(detuning, rabi, phase)
    @match (detuning, rabi, phase) begin
        (nothing, nothing, _) => error("at least one of rabi frequency or detuning parameter is required")
        (nothing, _, _) => (nothing, rabi, phase)
        (_, nothing, nothing) => (detuning, nothing, nothing)
        # if phase is specified, we set rabi to to its default value
        (_, nothing, _) => (detuning, 1.0, phase)
        _ => (detuning, rabi, phase)
    end
end

@option "rydberg" struct RydbergPulse
    detuning::PulseParameter = nothing
    rabi::PulseParameter = nothing
    phase::PulseParameter = nothing

    function RydbergPulse(detuning, rabi, phase)
        new(check_pulse_parameters(detuning, rabi, phase)...)
    end
end

@option "hyperfine" struct HyperfinePulse
    detuning::PulseParameter = nothing
    rabi::PulseParameter = nothing
    phase::PulseParameter = nothing

    function HyperfinePulse(detuning, rabi, phase)
        new(check_pulse_parameters(detuning, rabi, phase)...)
    end
end

@option struct Pulse
    atom_indices::Maybe{Vector{Int}} = nothing
    duration::Float64
    rydberg::Maybe{RydbergPulse} = nothing
    hyperfine::Maybe{HyperfinePulse} = nothing

    function Pulse(atom_indices, duration, rydberg, hyperfine)
        rydberg === hyperfine === nothing && throw(ValidationError("must specify rydberg or hyperfine"))
        new(atom_indices, duration, rydberg, hyperfine)
    end
end

@option struct PulseJob
    backend_name::String = "emulator"
    version::VersionNumber = v"0.1.0"
    job_id::UUID = uuid1()
    nshots::Int = 1000
    lattice_angle::Float64 = 90
    lattice_constant_a::Float64 = 3 * 0.844
    lattice_constant_b::Float64 = 3 * 0.844 * sin(lattice_angle)
    positions::Maybe{Vector{Tuple{Int, Int}}} = nothing
    natoms::Maybe{Int} = nothing
    pulses::Vector{Pulse}

    function PulseJob(backend_name, version, job_id, nshots,
            lattice_angle, lattice_constant_a, lattice_constant_b,
            positions, natoms, pulses)
        d = 0.844
        @assert 0 ≤ lattice_angle ≤ 180
        @assert lattice_constant_a ≥ 3d
        @assert lattice_constant_b * sind(lattice_angle) ≥ 3d
        natoms === nothing ||  @assert natoms ≤ 128
        # TODO: check lattice_constant_b is multiple of 0.1d
        new(backend_name, version, job_id, nshots,
            lattice_angle, lattice_constant_a, lattice_constant_b,
            positions, natoms, pulses)
    end
end
