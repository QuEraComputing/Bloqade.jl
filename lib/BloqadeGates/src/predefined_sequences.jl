function local_CkZ(atoms, ctrls, locs)
    n = length(atoms)
    @assert ctrls ⊆ 1:n
    @assert locs ⊆ 1:n
    @assert isempty(locs ∩ ctrls)
    seq = [local_single_qubit_gate(atoms, [ctrls; locs], X)]
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

function local_CkNOT(atoms, ctrls, locs)
    seq = [local_single_qubit_gate(atoms, locs, H)]
    append!(seq, local_CkZ(atoms, ctrls, locs))
    push!(seq, local_single_qubit_gate(atoms, locs, H))
    return seq
end
