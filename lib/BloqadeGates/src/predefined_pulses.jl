function local_single_qubit_gate(atoms, locs::Vector{<:Integer}, gate::AbstractBlock{2}; backend = KrylovEvolution, step = 1e-2)
    n = length(atoms)
    mask = zeros(n)
    mask[locs] .= 1
    Ω_hf, ϕ_hf, Δ_hf, t = single_qubit_gate_params(gate)
    return local_pulse(atoms, mask, t; Ω_hf, ϕ_hf, Δ_hf, backend, step)
end

function global_single_qubit_gate(atoms, gate::AbstractBlock{2}; backend = KrylovEvolution, step = 1e-2)
    Ω_hf, ϕ_hf, Δ_hf, t = single_qubit_gate_params(gate)
    return global_pulse(atoms, t; Ω_hf, ϕ_hf, Δ_hf, backend, step)
end

function single_qubit_gate_params(gate::AbstractBlock{2})
    @assert nqubits(gate) == 1
    @assert isunitary(gate)

    # Trivial Gates
    (gate isa Union{I2Gate, PhaseGate, TrivialGate}) && return (1.0, 0.0, 0.0, 0.0)

    # Constant Gates
    (gate isa XGate) && return (1.0, 0.0, 0.0, π)
    (gate isa YGate) && return (1.0, -π/2, 0.0, π)
    (gate isa ZGate) && return (0.0, 0.0, 1.0, π)
    (gate isa HGate) && return (sqrt(1/2), 0.0, sqrt(1/2), π)   # H = 1/√2 * (X + Z)
    (gate isa ConstGate.SGate) && return (0.0, 0.0, 1.0, π/2)
    (gate isa TGate) && return (0.0, 0.0, 1.0, π/4)
    (gate isa ConstGate.SdagGate) && return (0.0, 0.0, 1.0, 3π/2)
    (gate isa ConstGate.TdagGate) && return (0.0, 0.0, 1.0, 7π/4)

    # Rotation Gates
    (gate isa ShiftGate) && return (0.0, 0.0, 1, rem2pi(gate.theta, RoundDown))
    if gate isa RotationGate
        axis = gate.block
        θ = gate.theta
        (axis isa XGate) && return (1.0, 0.0, 0.0, rem2pi(θ, RoundDown))
        (axis isa YGate) && return (1.0, -π/2, 0.0, rem2pi(θ, RoundDown))
        (axis isa ZGate) && return (0.0, 0.0, 1.0, rem2pi(θ, RoundDown))
    end

    # Arbitrary 1-qubit Gates
    U = Matrix(gate)
    Ht = log(U) / -im
    n = real.([1/2*tr(Ht*mat(σ)) for σ in [X, Y, Z]])
    Ωt = 2*sqrt(n[1]^2 + n[2]^2)
    ϕ = -atan(n[2], n[1])
    Δt = 2*n[3]
    t = isapprox(Ωt, 0; atol = 1e-8) ? 1 : Ωt
    Ω, Δ = [Ωt, Δt]/t   # Ω = 1 if non-zero
    return Ω, ϕ, Δ, t
end
