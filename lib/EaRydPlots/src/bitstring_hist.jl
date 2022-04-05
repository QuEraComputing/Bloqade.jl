function find_largest(probs, num)
    if length(probs) ≤ num
        return 1:length(probs)
    end

    perm = sortperm(probs, rev=true)
    indices = perm[1:num]
    return indices
end

function bitstring_hist(r::ArrayReg; nlargest::Int, title="", kw...)
    ps = probs(r)
    indices = find_largest(ps, nlargest)

    obj = plt.bar(
        1:length(indices), ps[indices],
        kw...
    )
    plt.xlabel("bitstring")
    plt.ylabel("probability")
    plt.xticks(
        1:length(indices),
        string.(indices.-1;base=2, pad=nqubits(r)),
        rotation=60,
    )
    return obj
end

function bitstring_hist(r::SubspaceArrayReg; nlargest::Int, title="", kw...)
    state = statevec(r)
    probs = abs2.(state)
    indices = find_largest(probs, nlargest)

    obj = plt.bar(1:length(indices), probs[indices]; kw...)
    plt.xlabel("bitstring")
    plt.ylabel("probability")
    plt.xticks(
        1:length(indices),
        string.(space(r).subspace_v[indices];base=2, pad=nqubits(r)),
        rotation=60,
    )
    return obj
end
