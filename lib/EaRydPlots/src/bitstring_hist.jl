function find_largest(probs, num)
    if length(probs) ≤ num
        return 1:length(probs)
    end

    perm = sortperm(probs, rev=true)
    indices = perm[1:num]
    return indices
end

function bitstring_hist!(ax, r::ArrayReg; nlargest::Int, title="", kw...)
    ps = probs(r)
    indices = find_largest(ps, nlargest)

    obj = ax.bar(
        1:length(indices), ps[indices],
        kw...
    )
    ax.set_xlabel("bitstring")
    ax.set_ylabel("probability")
    ax.set_xticks(
        1:length(indices),
        string.(indices.-1;base=2, pad=nqubits(r)),
        rotation=60,
    )
    return obj
end

function bitstring_hist!(ax, r::SubspaceArrayReg; nlargest::Int, title="", kw...)
    state = statevec(r)
    probs = abs2.(state)
    indices = find_largest(probs, nlargest)

    obj = ax.bar(1:length(indices), probs[indices]; kw...)
    ax.set_xlabel("bitstring")
    ax.set_ylabel("probability")
    ax.set_xticks(
        1:length(indices),
        string.(space(r).subspace_v[indices];base=2, pad=nqubits(r)),
        rotation=60,
    )
    return obj
end

function bitstring_hist(r; kw...)
    fig, ax = plt.subplots()
    bitstring_hist!(ax, r; kw...)
    return fig
end
