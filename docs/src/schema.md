# Schema

Bloqade contains its own schema used to represent Hamiltonians in an IR (Intermediate Representation) that can then be executed via simulator or converted to and from other formats. Furthermore, tools such as [`hardware_transform`](@ref) and [`validate`](@ref) are available to check that user-defined Hamiltonians are capable of being executed on hardware and if not, transform them to be able to do so.

## Transforming Hamiltonians

We start with creating a hamiltonian:

```@repl transform_validation_example
using Bloqade, BloqadeSchema
Δ = constant(;duration=1.1, value=1.2*2π); # create detuning waveform
Ω = linear_ramp(duration=1.1, start_value=0.0, stop_value=2π*1.0); # create Rabi frequency waveform
ϕ = Waveform(t->2.2*2π*sin(t)^2, duration=1.1); # create phase waveform
atoms = generate_sites(ChainLattice(), 4, scale=1.0); # provide lattice geometry
h = rydberg_h(atoms; Δ = Δ, Ω = Ω, ϕ = ϕ) # put it all together to create a hamiltonian
```

To transform the Hamiltonian into something the hardware is capable of supporting, we can pass it through [`hardware_transform`](@ref). `hardware_transform` accepts information from [`get_device_capabilities`](@ref) (already called as a default argument) which provides information on the machine's capabilities and returns the transformed hamiltonian along with additional information  regarding the difference between the originally defined lattice geometry and waveforms versus their transformed versions.

```@repl transform_validation_example
transformed_h, transform_info = hardware_transform(h);
transformed_h
transform_info
dump(transform_info) # the returned HardwareTransformInfo has fields containing corresponding transformation information
```

## Validating Hamiltonians


We can see if this hamiltonian or any other hamiltonian we create will run on hardware with the help of [`validate`](@ref). 
```@repl transform_validation_example
validate(transformed_h)
```

## Reference

```@docs
get_device_capabilities
get_device_capabilities_SI
get_rydberg_capabilities
hardware_transform_Ω
hardware_transform_ϕ
hardware_transform_Δ
hardware_transform_atoms
hardware_transform
validate
to_json
to_json_no_validation
to_dict
to_schema
to_schema_no_validation
from_json
from_dict
from_schema
execute
TaskSpecification
TaskOutput
ValidationViolations
```