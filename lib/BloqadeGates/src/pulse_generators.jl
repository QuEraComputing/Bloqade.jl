function apply_mask(p, mask::Vector{<:Real})
    p isa Real && return p * mask
    p isa Function && return [(x->p(x)*mask[i]) for i in eachindex(mask)]
    error("Pulse parameter should be a real number or a function.")
end

function local_hamiltonian(atoms, mask::Vector{<:Real}; 
        Ω_hf = nothing, ϕ_hf = nothing, Δ_hf = nothing, 
        Ω_r = nothing, ϕ_r = nothing, Δ_r = nothing)
    !isnothing(Ω_hf) && (Ω_hf = apply_mask(Ω_hf, mask))
    !isnothing(ϕ_hf) && (ϕ_hf = apply_mask(ϕ_hf, mask))
    !isnothing(Δ_hf) && (Δ_hf = apply_mask(Δ_hf, mask))
    !isnothing(Ω_r) && (Ω_r = apply_mask(Ω_r, mask))
    !isnothing(ϕ_r) && (ϕ_r = apply_mask(ϕ_r, mask))
    !isnothing(Δ_r) && (Δ_r = apply_mask(Δ_r, mask))
    return rydberg_h_3(atoms; Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r)
end

function local_pulse(atoms, mask::Vector{<:Real}, t::Real; backend = KrylovEvolution, step = 1e-2,
        Ω_hf = nothing, ϕ_hf = nothing, Δ_hf = nothing, 
        Ω_r = nothing, ϕ_r = nothing, Δ_r = nothing)
    h = local_hamiltonian(atoms, mask; Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r)
    return RydbergPulse(h, t; backend, step)
end

function global_hamiltonian(atoms; 
        Ω_hf = nothing, ϕ_hf = nothing, Δ_hf = nothing, 
        Ω_r = nothing, ϕ_r = nothing, Δ_r = nothing)
    return rydberg_h_3(atoms; Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r)
end

function global_pulse(atoms, t::Real; backend = KrylovEvolution, step = 1e-2,
        Ω_hf = nothing, ϕ_hf = nothing, Δ_hf = nothing, 
        Ω_r = nothing, ϕ_r = nothing, Δ_r = nothing)
    h = global_hamiltonian(atoms; Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r)
    return RydbergPulse(h, t; backend, step)
end

function single_site_pulse(atoms, j::Integer, t::Real; backend = KrylovEvolution, step = 1e-2,
        Ω_hf = nothing, ϕ_hf = nothing, Δ_hf = nothing, 
        Ω_r = nothing, ϕ_r = nothing, Δ_r = nothing)
    n = length(atoms)
    mask = zeros(n)
    mask[j] = 1
    return local_pulse(atoms, mask, t; Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r, backend, step)
end