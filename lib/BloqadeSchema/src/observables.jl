"""
    rydberg_density(task_res::AnalogHamiltonianSimulationQuantumTaskResult) -> Real

Return the rydberg density at each site given an `AnalogHamiltonianSimulationTaskResult`.
"""
rydberg_density(task_res::Braket.AnalogHamiltonianSimulationQuantumTaskResult) = sum(map(x -> 1 .- x.post_sequence, task_res.measurements)) / length(task_res.measurements)

"""
    rydberg_density(task_res::AnalogHamiltonianSimulationQuantumTaskResult, i::Int) -> Real


Return the rydberg density at site `i` given an `AnalogHamiltonianSimulationTaskResult`.
"""
rydberg_density(task_res::Braket.AnalogHamiltonianSimulationQuantumTaskResult, i::Int) = sum(map(x -> 1 .- x.post_sequence[i], task_res.measurements)) / length(task_res.measurements)