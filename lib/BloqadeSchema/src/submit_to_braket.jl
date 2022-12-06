using BloqadeSchema
using BloqadeExpr
using AWS

"""
    submit_to_braket(h::BloqadeExpr.Hamiltonian, n_shots::Int; <keyword arguments>)

Submits a `BloqadeExpr.Hamiltonian` instance to Braket with `n_shots` defining the number of times
the Hamiltonian should be executed. 



# Keyword Arguments
- `arn="arn:aws:braket:us-east-1::device/qpu/quera/Aquila"`: ARN for the machine
- `region="us-east-1"`: AWS Region machine is located in
- `credentials::Union{AWSCredentials, Nothing}=nothing`: `AWS.AWSCredentials` instance you can create to login.

Credentials can be passed in explicitly through an `AWS.AWSCredentials` struct or by passing in 
`nothing`, in which case credentials will be sought after in the shell environment variables. 

# Keyword Arguments
- `arn="arn:aws:braket:us-east-1::device/qpu/quera/Aquila"`: ARN for the machine
- `region="us-east-1"`: AWS Region machine is located in
- `credentials::Union{AWSCredentials, Nothing}=nothing`: `AWS.AWSCredentials` instance you can create to login.
"""
function submit_to_braket(h::BloqadeExpr.Hamiltonian, 
                          n_shots::Int, 
                          device_capabilities = BloqadeSchema.get_device_capabilities();
                            arn="arn:aws:braket:us-east-1::device/qpu/quera/Aquila",
                            region="us-east-1",
                            credentials::Union{AWSCredentials, Nothing}=nothing
                        ) 
    # is a BloqadeSchema.TaskSpecification instance, can call the other function
    translation_params = BloqadeSchema.SchemaTranslationParams(
        n_shots,
        device_capabilities
    )
    return submit_to_braket(h, translation_params; arn=arn, region=region, credentials=credentials)
end

"""
    submit_to_braket(h::BloqadeExpr.Hamiltonian, translation_params::BloqadeSchema.SchemaTranslationParams; <keyword arguments>)

Submits a `BloqadeExpr.Hamiltonian` instance to Braket with `BloqadeSchema.SchemaTranslationParams` containing the number of shots and 
device capabilities, returning an `AWS.AwsQuantumTask` and `Bloqade.HardwareTransformInfo` upon
converting the Hamiltonian to one the hardware can execute.

Credentials can be passed in explicitly through an `AWS.AWSCredentials` struct or by passing in 
`nothing`, in which case credentials will be sought after in the shell environment variables. 

# Keyword Arguments
- `arn="arn:aws:braket:us-east-1::device/qpu/quera/Aquila"`: ARN for the machine
- `region="us-east-1"`: AWS Region machine is located in
- `credentials::Union{AWSCredentials, Nothing}=nothing`: `AWS.AWSCredentials` instance you can create to login.
"""
function submit_to_braket(h::BloqadeExpr.Hamiltonian,
                          translation_params::BloqadeSchema.SchemaTranslationParams;
                            arn="arn:aws:braket:us-east-1::device/qpu/quera/Aquila",
                            region="us-east-1",
                            credentials::Union{AWSCredentials, Nothing}=nothing
                        )

    h_transformed, transform_info = hardware_transform(h)
    bloqade_ir = BloqadeSchema.to_schema(h_transformed, translation_params)
    task = submit_to_braket(bloqade_ir; arn=arn, region=region, credentials=credentials)
    return task, transform_info
end

"""
    sumbit_to_braket(ts:BloqadeSchema.TaskSpecification; <keyword arguments>)

Submits a `BloqadeSchema.TaskSpecification` instance to Braket, returning an `AWS.AwsQuantumTask`.

Credentials can be passed in explicitly through an `AWS.AWSCredentials` struct or by passing in 
`nothing`, in which case credentials will be sought after in the shell environment variables. 

# Keyword Arguments
- `arn="arn:aws:braket:us-east-1::device/qpu/quera/Aquila"`: ARN for the machine
- `region="us-east-1"`: AWS Region machine is located in
- `credentials::Union{AWSCredentials, Nothing}=nothing`: `AWS.AWSCredentials` instance you can create to login.
"""
function submit_to_braket(ts::BloqadeSchema.TaskSpecification; 
        arn="arn:aws:braket:us-east-1::device/qpu/quera/Aquila", 
        region="us-east-1",
        credentials::Union{AWSCredentials, Nothing}=nothing)
    # may need region provided, don't think we have to worry about output

    aws_config = if !isnothing(credentials)
        AWSConfig(;creds=credentials, region=region)
    else
        # default to environment variables
        AWSConfig(;region=region)
    end
    
    qpu = Braket.AwsDevice(arn; config = aws_config)
    braket_ir = to_braket_ahs_ir(ts)
    task = qpu(braket_ir;shots=ts.nshots)
    
    return task
end