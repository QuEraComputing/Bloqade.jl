# Schema

Bloqade contains its own schema used to represent Hamiltonians in an IR (Intermediate Representation) that can then be executed via simulator or converted to and from other formats. Furthermore, tools such as [`hardware_transform`](@ref) and [`validate`](@ref) are available to check that user-defined Hamiltonians are capable of being executed on hardware and if not, transform them to be able to do so.

!!! warning "3-Level Support"
    The schema and conversion functionalities are currently not available for 3-level Hamiltonians

## Transforming Hamiltonians

We start with creating a hamiltonian:

```@repl schema_example
using Bloqade, BloqadeSchema
Δ = constant(;duration=1.1, value=1.2*2π); 
Ω = linear_ramp(duration=1.1, start_value=0.0, stop_value=2π*1.0); 
ϕ = Waveform(t->2.2*2π*sin(t)^2, duration=1.1);
atoms = generate_sites(ChainLattice(), 4, scale=1.0);
h = rydberg_h(atoms; Δ = Δ, Ω = Ω, ϕ = ϕ)
```

To transform the Hamiltonian into something the hardware is capable of supporting, we can pass it through [`hardware_transform`](@ref).

`hardware_transform` accepts information from [`get_device_capabilities`](@ref) (already called as a default argument) which provides information on the machine's capabilities and returns the transformed hamiltonian along with additional information regarding the difference (error) between the originally defined lattice geometry and waveforms versus their transformed versions through a [`HardwareTransformInfo`](@ref) instance.


```@repl schema_example
transformed_h, transform_info = hardware_transform(h);
transformed_h
transform_info
dump(transform_info)
```

## Validating Hamiltonians


We can see if this hamiltonian or any other hamiltonian we create will run on hardware with the help of [`validate`](@ref). 
```@repl schema_example
validate(transformed_h)
```

In this case, the waveforms have been successfully transformed but there are still some issues with the atom positions. We can rescale their positions, regenerate the Hamiltonian and validate again to make sure the changes are correct.

```@repl schema_example
atoms = generate_sites(ChainLattice(), 4, scale=4.0);
fixed_h = rydberg_h(atoms; Δ = Δ, Ω = Ω, ϕ = ϕ) # Keep older waveforms with new atom geometry
transformed_h, _ = hardware_transform(fixed_h)
validate(transformed_h)
```

No violations are present meaning the new Hamiltonian has passed validation and can now be executed on hardware.

## Converting Between Formats

You can convert the Hamiltonian to and from:
* JSON Object format
* Julia Dictionary representation
* Schema representation
to store the Hamiltonian for other applications.

!!! note "Internal Validation"
    By default, all conversion functions invoke [`validate`](@ref) internally to ensure the Hamiltonian is capable of being run on hardware. If a violation is detected, the Hamiltonian is not converted. For JSON and Schema representations, this can be bypassed by invoking the "no validation" variants [`to_schema_no_validation`](@ref) and [`to_json_no_validation`](@ref) respectively. This bypass ability is not available for dictionary representation.

### Schema

[`to_schema`](@ref) allows you to convert a Hamiltonian to Bloqade's native schema format ([`TaskSpecification`](@ref)) along with an optional
[`SchemaTranslationParams`](@ref) argument to specify the number of shots and device capabilities for validation.

```@repl schema_example
h_schema = to_schema(transformed_h)
from_schema(h_schema) # to convert back to Hamiltonian
```

### JSON
[`to_json`](@ref) allows you to convert a Hamiltonian to a JSON Object along with an optional
[`SchemaTranslationParams`](@ref) argument to specify the number of shots and device capabilities for validation.
```@repl schema_example
h_json = to_json(transformed_h)
from_json(h_json) # to convert back to Hamiltonian
```

### Dictionary

[`to_dict`](@ref) allows you to convert a Hamiltonian to a Julia dictionary along with an optional 
[`SchemaTranslationParams`](@ref) argument to specify the number of shots and device capabilities for validation.
```@repl schema_example
h_dict = to_dict(transformed_h)
from_dict(h_dict) # to convert back to Hamiltonian
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
HardwareTransformInfo
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