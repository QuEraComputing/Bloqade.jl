function local_CkZ(atoms, ctrls::Vector{<:Integer}, locs::Vector{<:Integer}; 
        backend = KrylovEvolution, step = 1e-2)
    n = length(atoms)
    @assert !isempty(ctrls) && !isempty(locs)
    @assert ctrls ⊆ 1:n
    @assert locs ⊆ 1:n
    @assert isempty(locs ∩ ctrls)
    seq = chain(n, local_single_qubit_gate(atoms, [ctrls; locs], X; backend, step))
    append!(seq, [single_site_pulse(atoms, ctrl, π; Ω_r = 1.0, ϕ_r = 0.0, Δ_r = 0.0, backend, step) for ctrl in ctrls])
    for loc in locs
        mask = zeros(n)
        mask[loc] = 1
        push!(seq, local_pulse(atoms, mask, 2π; Ω_r = 1.0, ϕ_r = 0.0, Δ_r = 0.0, backend, step))
    end
    append!(seq, reverse!([single_site_pulse(atoms, ctrl, π; Ω_r = 1.0, ϕ_r = 0.0, Δ_r = 0.0, backend, step) for ctrl in ctrls]))
    push!(seq, local_single_qubit_gate(atoms, [ctrls; locs], X; backend, step))
    return seq
end

function local_CkNOT(atoms, ctrls::Vector{<:Integer}, locs::Vector{<:Integer}; 
        backend = KrylovEvolution, step = 1e-2)
    seq = chain(length(atoms), local_single_qubit_gate(atoms, locs, H; backend, step))
    append!(seq, local_CkZ(atoms, ctrls, locs; backend, step))
    push!(seq, local_single_qubit_gate(atoms, locs, H; backend, step))
    return seq
end

function global_levine_pichler(atoms; Ω_r = 1.0, backend = KrylovEvolution, step = 1e-2)
    ϕ_r = 3.90242
    Δ_r = 0.377371*Ω_r
    τ = 4.29268/Ω_r
    seq = chain(
        global_pulse(atoms, τ; Ω_r, ϕ_r = 0.0, Δ_r, backend, step),
        global_pulse(atoms, τ; Ω_r, ϕ_r, Δ_r, backend, step)
    )
    return seq
end
