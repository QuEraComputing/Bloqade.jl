function find_largest(probs, num)
    if length(probs) ≤ num
        return 1:length(probs)
    end

    perm = sortperm(probs, rev = true)
    indices = perm[1:num]
    return indices
end

"""
    bitstring_hist!(ax, register; nlargest::Int, title="", kw...)

Plot the bitstring histogram.

# Arguments

- `ax`: the axis object from `CairoMakie`.
- `register`: the register to plot.

# Keyword Arguments

- `nlargest`: plot the first `nlargest` bitstrings.
- `title`: title of the plot.
- `kw`: other keyword arguments supported by `CairoMakie.bars!` function.
"""
function bitstring_hist!(ax, r::ArrayReg; nlargest::Int, title = "", kw...)
    ps = probs(r)
    indices = find_largest(ps, nlargest)

    obj = barplot!(ax, 1:length(indices), ps[indices]; kw...)
    return obj
end

function bitstring_hist!(ax, r::SubspaceArrayReg; nlargest::Int, title = "", kw...)
    state = statevec(r)
    probs = abs2.(state)
    indices = find_largest(probs, nlargest)

    obj = barplot!(ax, 1:length(indices), probs[indices]; kw...)
    return obj
end

"""
    bitstring_hist(r; kw...)

Plot the bitstring histogram.

# Arguments

- `register`: the register to plot.

# Keyword Arguments

- `nlargest`: plot the first `nlargest` bitstrings.
- `title`: title of the plot.
- `kw...`: other keyword arguments supported by `CairoMakie.bars!` function.
"""
function bitstring_hist(r; nlargest::Int, kw...)
    fig = CairoMakie.Figure()
    xticks_string = _xticks_string(r, nlargest)
    ax = Axis(fig[1, 1], xlabel="bitstring", ylabel="probability", xticks = (1:length(xticks_string), xticks_string), xticklabelrotation=π/3)

    bitstring_hist!(ax, r; nlargest, kw...)
    return fig
end
function _xticks_string(r::ArrayReg, nlargest::Int)
    ps = probs(r)
    indices = find_largest(ps, nlargest)
    return string.(indices .- 1; base = 2, pad = nqubits(r))
end
function _xticks_string(r::SubspaceArrayReg, nlargest::Int)
    state = statevec(r)
    probs = abs2.(state)
    indices = find_largest(probs, nlargest)
    return string.(space(r).subspace_v[indices]; base = 2, pad = nqubits(r))
end