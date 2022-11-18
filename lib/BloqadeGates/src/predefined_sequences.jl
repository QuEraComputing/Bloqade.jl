function local_CkZ(atoms, ctrls::Vector{<:Integer}, locs::Vector{<:Integer})
    n = length(atoms)
    @assert !isempty(ctrls) && !isempty(locs)
    @assert ctrls ⊆ 1:n
    @assert locs ⊆ 1:n
    @assert isempty(locs ∩ ctrls)
    seq = chain(n, local_single_qubit_gate(atoms, [ctrls; locs], X))
    append!(seq, [single_site_pulse(atoms, ctrl, 1.0, 0.0, 0.0, π; pulse_type = :rydberg) for ctrl in ctrls])
    for loc in locs
        mask = zeros(n)
        mask[loc] = 1
        push!(seq, local_pulse(atoms, mask, 1.0, 0.0, 0.0, 2π; pulse_type = :rydberg))
    end
    append!(seq, reverse!([single_site_pulse(atoms, ctrl, 1.0, 0.0, 0.0, π; pulse_type = :rydberg) for ctrl in ctrls]))
    push!(seq, local_single_qubit_gate(atoms, [ctrls; locs], X))
    return seq
end

function local_CkNOT(atoms, ctrls::Vector{<:Integer}, locs::Vector{<:Integer})
    seq = chain(length(atoms), local_single_qubit_gate(atoms, locs, H))
    append!(seq, local_CkZ(atoms, ctrls, locs))
    push!(seq, local_single_qubit_gate(atoms, locs, H))
    return seq
end

function global_levine_pichler(atoms; Ω_r = 1.0)
    ϕ_r = 3.90242
    Δ_r = 0.377371*Ω_r
    τ = 4.29268/Ω_r
    seq = chain(
        global_pulse(atoms, Ω_r, 0.0, Δ_r, τ; pulse_type = :rydberg),
        global_pulse(atoms, Ω_r, ϕ_r, Δ_r, τ; pulse_type = :rydberg)
    )
    return seq
end
