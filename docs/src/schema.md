# [Interacting with Neutral Atom Hardware](@id schema)

Bloqade contains its own schema used to represent Hamiltonians in an IR (Intermediate Representation) that can then be executed via simulator/hardware as well as converted to and from other formats. Furthermore, tools such as [`hardware_transform`](@ref) and [`validate`](@ref) are available to check that user-defined Hamiltonians are capable of being executed on hardware and if not, transform them to be able to do so.

!!! warning "3-Level Support"
    The schema and conversion capabilities are currently not available for 3-level Hamiltonians

## Transforming Hamiltonians to Hardware Compatible Form

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

!!! warning "Limitations on Atom Position Transformations"
    While `hardware_transform` may attempt to adjust atom positions so that they conform to hardware position resolution capabilities, the function will NOT move atoms such that they satisfy minimum spacing constraints. The `validate` function presented later will explicitly indicate which atoms are in violation of the position constraints but will require the user to make the necessary changes.

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
* Amazon Braket representation
to store the Hamiltonian for other applications or execute it in Bloqade/Neutral Atom hardware. 

!!! note "Internal Validation"
    By default, all conversion functions invoke [`validate`](@ref) internally to ensure the Hamiltonian is capable of being run on hardware. If a violation is detected, the Hamiltonian is not converted. For JSON and Schema representations, this can be bypassed by invoking the "no validation" variants [`to_schema_no_validation`](@ref) and [`to_json_no_validation`](@ref) respectively. This bypass ability is not available for dictionary representation and Amazon Braket representation (which requires the Hamiltonian is already in [`TaskSpecification`](@ref) format).

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

### Amazon Braket
[`to_braket_ahs_ir`](@ref) allows you to convert a [`TaskSpecification`](@ref) instance into an [`Braket.IR.AHSProgram`](https://github.com/awslabs/Braket.jl/blob/main/src/raw_schema.jl#L604) that can be submitted for execution on Neutral Atom hardware such as [QuEra's 256-qubit Aquila machine](https://www.quera.com/aquila) via [Amazon Braket](https://aws.amazon.com/braket/).
```@repl schema_example
h_schema = to_schema(transformed_h)
h_braket = to_braket_ahs_ir(h_schema)
```

!!! info "One-Way Conversion"
    Unlike the previous conversion functions, `from_braket_ahs_ir` does not exist meaning it is not possible to obtain a [`TaskSpecification`](@ref) from a `Braket.IR.AHSProgram`.

!!! warning "No Validation"
    [`to_braket_ahs_ir`](@ref) does no validation on its input. It is assumed that the [`TaskSpecification`](@ref) given to it is already validated. Therefore, it is possible to produce a `Braket.IR.AHSProgram` that is incapable of being run on Neutral Atom hardware and may cause the Braket API to reject it.

## Submitting to Amazon Braket

To submit to Neutral Atom hardware on Amazon Braket, Bloqade provides [`submit_to_braket`](@ref) which can submit BOTH the native Bloqade representation of Hamiltonians (`BloqadeExpr.RydbergHamiltonian`) as well as the [`TaskSpecification`](@ref) representation.

!!! warning "Implicit Transformation and No Validation for TaskSpecification"
    For any `BloqadeExpr.RydbergHamiltonian` passed in, [`hardware_transform`](@ref) is invoked to ensure it is compatible with hardware. On the other hand, [`TaskSpecification`](@ref) types are assumed to already be valid.


[`submit_to_braket`](@ref) requires that AWS credentials are given either explicitly through an [`AWS.AWSCredentials`](https://github.com/JuliaCloud/AWS.jl/blob/master/src/AWSCredentials.jl#L32) type or by setting the environment variables in the shell running Bloqade with the credentials. The credentials (and instructions for setting environment variables!) can be found through your AWS account's "Command line or programmatic access" option.

Let us try to submit the Hamiltonian we made earlier. We remind ourselves that our Hamiltonian is currently the following:

```@repl schema_example
fixed_h
```

(Recall we had to modify the atom positions in order to pass validation!)

Now we define the number of shots (how many times the Hamiltonian will be executed on hardware) as well as the credentials, allowing [`submit_to_braket`](@ref) to automatically handle transforming the Hamiltonian to fit within hardware capabilities. By default, [`submit_to_braket`](@ref) will submit to QuEra's Aquila Neutral Atom hardware and take into account its capabilities for Hamiltonian transformation.

```julia
using AWS
access_key_id = "your_access_key_id"
secret_key = "your_secret_key"
token = "your_token"
credentials = AWS.AWSCredentials(access_key_id, secret_key, token)
task = submit_to_braket(fixed_h, 100; credentials=credentials)
```

If submission was successful you will see something like 
```julia
AwsQuantumTask(...)
```
in the REPL with the Task ARN (Amazon Resource Name) as a string in the parentheses.

## Inspecting Results from Braket

To see the status of our task we can use `state` from the `Braket.jl` package.

```julia
using Braket
Braket.state(task)
```

`state` can return a `String` that is either: `"CANCELLED"`, `"FAILED"`, `"COMPLETED"`, `"QUEUED"`, or `"RUNNING"`.

To obtain results, the `result` function from `Braket.jl` can be used

```julia
result = Braket.result(task)
```

!!! info "Braket.result is Blocking"
    Per the docstring for [`result`](https://github.com/awslabs/Braket.jl/blob/main/src/task.jl#L292), the function is **BLOCKING** "until a result
    is available, in which case the result is returned, or the task enters a
    terminal state without a result (`"FAILED"` or `"CANCELLED"`)...".

To obtain the raw measurements (pre- and post-Hamiltonian application) of the atoms, the `get_measurements` function in `Braket.jl` can be used.

```julia
Braket.get_measurements(result)
```

## Reference

```@docs
get_device_capabilities
get_device_capabilities_SI
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
to_braket_ahs_ir
submit_to_braket
from_json
from_dict
from_schema
execute
TaskSpecification
TaskOutput
ValidationViolations
```
