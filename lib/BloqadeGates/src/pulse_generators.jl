function local_hamiltonian(atoms, mask, Ω, ϕ, Δ; pulse_type = :rydberg)
    Ω = Ω * mask
    ϕ = ϕ * mask
    Δ = Δ * mask
    pulse_type === :rydberg && (h = rydberg_h_3(atoms; Ω_r = Ω, ϕ_r = ϕ, Δ_r = Δ))
    pulse_type === :hyperfine && (h = rydberg_h_3(atoms; Ω_hf = Ω, ϕ_hf = ϕ, Δ_hf = Δ))
    return h
end

function local_pulse(atoms, mask, Ω, ϕ, Δ, t; 
    pulse_type = :rydberg, backend = KrylovEvolution, step = 1e-2)
    h = local_hamiltonian(atoms, mask, Ω, ϕ, Δ; pulse_type)
    return RydbergPulse(h, t; backend, step)
end

function global_hamiltonian(atoms, Ω, ϕ, Δ; pulse_type = :rydberg)
    pulse_type === :rydberg && (h = rydberg_h_3(atoms; Ω_r = Ω, ϕ_r = ϕ, Δ_r = Δ))
    pulse_type === :hyperfine && (h = rydberg_h_3(atoms; Ω_hf = Ω, ϕ_hf = ϕ, Δ_hf = Δ))
    return h
end

function global_pulse(atoms, Ω, ϕ, Δ, t; 
        pulse_type = :rydberg, backend = KrylovEvolution, step = 1e-2)
    h = global_hamiltonian(atoms, Ω, ϕ, Δ; pulse_type)
    return RydbergPulse(h, t; backend, step)
end
