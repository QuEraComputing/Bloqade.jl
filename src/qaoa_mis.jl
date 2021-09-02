export mean_rydberg, count_vertices, mean, gibbs_loss, logsumexp

"""
    count_vertices(config::Integer)

counter the number of vertices in a spin configuration.
"""
count_vertices(config::Integer) = count_ones(config)
count_vertices(config::AbstractVector) = count(isone, config)

"""
    mean_rydberg([f], reg_or_samples)

Mean size of vertex set.

# Arguments

- `f`: optional, postprocessing callback function `f(config) -> config`.
    The input `config` is an integer of type `Int`, the output
    `config` can be a type supports [`count_vertices`](@ref)
    e.g, an `AbstractVector` or an `Integer`.
- `reg_or_samples` can be a register (`Yao.ArrayReg` or [`RydbergReg`](@ref))
    or a list of measurement result (config) in `AbstractVector`.

# Example

To implement the postprocessing protocal in MIS experiment:

1. calculating `mean_rydberg` by first reducing the configuration
to independent set using [`to_independent_set`](@ref)
2. randomly adding vertices then pick the largest [`count_vertices`](@ref)
using [`add_random_vertices!`](@ref).

```julia
mean_rydberg(r) do config
    config = to_independent_set(graph, config)
    add_random_vertices!(config, graph, 10)
    return config
end
```

Or one can also just add vertice by atom order

```julia
mean_rydberg(r) do config
    config = to_independent_set(graph, config)
    add_vertices!(config, graph)
    return config
end
```
"""
function mean_rydberg end

mean_rydberg(x) = mean_rydberg(identity, x)

_config_amplitude(r::RydbergReg) = zip(vec(r.subspace), vec(r.state))

function _config_amplitude(r::Yao.ArrayReg)
    state = Yao.relaxedvec(r)
    return zip(0:length(state)-1, state)
end

function mean_rydberg(f, reg::Yao.AbstractRegister)
    mean_ryd = zero(real(eltype(reg.state)))
    for (c, amp) in _config_amplitude(reg)
        nvertices = count_vertices(f(c))
        mean_ryd += abs2(amp) * nvertices
    end
    return mean_ryd
end

function mean_rydberg(f, samples::AbstractVector)
    return mean(count_vertices ∘ f, samples)
end

"""
    gibbs_loss([f], reg_or_samples, α::Real)

The Gibbs loss for maximum independent set defined as

```math
L = -1/α \\log(\\langle ψ|\\exp(α \\sum(n))|ψ\\rangle),
```

where `n` is the vertex set size.

# Arguments

- `f`: optional, postprocessing callback function `f(config) -> config`.
    The input `config` is an integer of type `Int`, the output
    `config` can be a type supports [`count_vertices`](@ref)
    e.g, an `AbstractVector` or an `Integer`.
- `reg_or_samples` can be a register (`Yao.ArrayReg` or [`RydbergReg`](@ref))
    or a list of measurement result (config) in `AbstractVector`.
- `α::Real`: the parameter of Gibbs loss.
"""
gibbs_loss(reg_or_samples, α::Real) = gibbs_loss(identity, reg_or_samples, α)

function gibbs_loss(f, reg::Yao.AbstractRegister, α::Real)
    expected = sum(_config_amplitude(reg)) do (config, amp)
        abs2(amp) * exp(α * count_vertices(config))
    end
    return -log(expected) / α 
end

function logsumexp(x::AbstractArray)
    xmax = maximum(x)
    return log(sum(exp.(x .- xmax))) + xmax
end

function gibbs_loss(f, samples::AbstractVector, α::Real)
    expect = map(x->α * count_vertices(f(x)), samples)
    return -(logsumexp(expect) - log(length(samples)))/α
end

function gibbs_loss(f, samples::AbstractMatrix, α::Real)
    -(logsumexp(α .* sum(samples, dims=(2,))) - log(size(samples, 1)))/α
end

gibbs_loss(α::Real) = reg -> gibbs_loss(reg, α)
