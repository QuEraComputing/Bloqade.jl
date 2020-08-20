export qaoa_on_graph, mean_nv, count_vertices, mean, gibbs_loss, logsumexp

"""
    count_vertices(config::Integer)

counter the number of vertices in a spin configuration.
"""
count_vertices(config::Integer) = count_ones(config)

"""
    qaoa_on_graph(graph, ts, ϕs[, cache, state])

Execute qaoa circuit on a graph model with `simple_rydberg` hamiltonians
and return a [`RydbergReg`](@ref) register. One can preallocate
the cache and register state by feeding two more arugments. The state
has to be an `Array` or `CuArray`, feeding other types will cause undefined
behaviour.
"""
function qaoa_on_graph end

function qaoa_on_graph(graph, ts::AbstractVector{<:Real}, ϕs::AbstractVector, cache::EmulatorCache, state)
    hs = simple_rydberg.(nv(graph), ϕs)
    subspace = Subspace(graph)
    fill!(state, zero(eltype(ts)))
    state[1, :] .= 1
    reg = RydbergReg{nv(graph)}(state, subspace)
    # evolve
    emulate!(reg, ts, hs, cache)
    return reg
end

function qaoa_on_graph(graph, ts::AbstractVector{<:Real}, ϕs::AbstractVector)
    hs = simple_rydberg.(nv(graph), ϕs)
    subspace = Subspace(graph)
    cache = EmulatorCache(ts, hs, subspace)
    r = zero_state(Complex{eltype(ts)}, nv(graph), subspace)
    emulate!(r, ts, hs, cache)
    return r
end

"""
    mean_nv(reg)

Mean size of vertex set.
`reg` can be a measurement result or a register.
"""
function mean_nv(reg::RydbergReg)
    sum(t -> abs2(t[2]) * count_vertices(t[1]), zip(vec(reg.subspace), Yao.relaxedvec(reg)))
end

function mean_nv(reg::Yao.ArrayReg)
    return sum(enumerate(Yao.relaxedvec(reg))) do (c, amp)
        abs2(amp) * count_vertices(c - 1)
    end
end

mean_nv(samples::AbstractVector{<:BitStr}) = mean(count_vertices, samples)

"""
    gibbs_loss([reg::RydbergReg], α::Real)

The Gibbs loss for maximum independent set defined as
```math
L = -1/α \\log(\\langle ψ|\\exp(α \\sum(n))|ψ\\rangle),
```
where `n` is the vertex set size.
"""
function gibbs_loss(reg::RydbergReg, α::Real)
    expected = sum(zip(vec(reg.subspace), Yao.relaxedvec(reg))) do (config, amp)
        abs2(amp) * exp(α * count_vertices(config))
    end
    -log(expected)/α
end

function gibbs_loss(reg::Yao.ArrayReg, α::Real)
    expected = sum(enumerate(Yao.relaxedvec(reg))) do (config, amp)
        abs2(amp) * exp(α * count_vertices(config - 1))
    end
    return -log(expected) / α
end

function logsumexp(x::AbstractArray)
    xmax = maximum(x)
    log(sum(exp.(x .- xmax))) + xmax
end

function gibbs_loss(samples::AbstractVector{<:BitStr}, α::Real)
    -(logsumexp(α .* count_vertices.(samples)) - log(length(samples)))/α
end

function gibbs_loss(samples::AbstractMatrix, α::Real)
    -(logsumexp(α .* sum(samples, dims=(2,))) - log(size(samples, 1)))/α
end

gibbs_loss(α::Real) = reg -> gibbs_loss(reg, α)
