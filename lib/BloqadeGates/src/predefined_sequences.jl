function local_CkZ(atoms, ctrls, locs)
    n = length(atoms)
    @assert ctrls ⊆ 1:n
    @assert locs ⊆ 1:n
    @assert isempty(locs ∩ ctrls)
    seq = [single_site_pi_pulse(atoms, i, :hyperfine) for i = [ctrls; locs]]
    append!(seq, [single_site_pi_pulse(atoms, ctrl, :rydberg) for ctrl in ctrls])
    for loc in locs
        mask = zeros(n)
        mask[loc] = 1
        push!(seq, local_pulse(atoms, mask, 1.0, 0.0, 0.0, 2π; pulse_type = :rydberg))
    end
    append!(seq, reverse!([single_site_pi_pulse(atoms, ctrl, :rydberg) for ctrl in ctrls]))
    append!(seq, [single_site_pi_pulse(atoms, i, :hyperfine) for i = [ctrls; locs]])
    return seq
end

function local_CkNOT(atoms, ctrls, locs)
    seq = [local_hadamard(atoms, loc) for loc in locs]
    append!(seq, local_CkZ(atoms, ctrls, locs))
    append!(seq, [local_hadamard(atoms, loc) for loc in locs])
    return seq
end
