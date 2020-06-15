export qaoa_on_graph, mean_nv, count_vertices, mean, soft_misloss, logsumexp

"""
    count_vertices(config::Integer)

counter the number of vertices in a spin configuration.
"""
count_vertices(config::Integer) = count_ones(config)

"""
    qaoa_on_graph(graph, ϕs::AbstractVector, ts::AbstractVector)

Execute qaoa circuit on a graph model and return a `RydbergReg` register.
"""
function qaoa_on_graph(graph, ϕs::AbstractVector, ts::AbstractVector)
    hs = SimpleRydberg.(ϕs)
    # prepair a zero state
    subspace_v = subspace(graph)
    st = zeros(ComplexF64, length(subspace_v)); st[1] = 1
    reg = RydbergReg{nv(graph)}(st, subspace_v)

    # evolve
    reg |> QAOA{nv(graph)}(subspace_v, hs, ts)
    return reg
end

"""
    mean_nv(reg)

Mean size of vertex set.
`reg` can be a measurement result or a register.
"""
function mean_nv(reg::RydbergReg)
    sum(t -> abs2(t[2]) * count_vertices(t[1]), zip(reg.subspace, relaxedvec(reg)))
end

mean_nv(samples::AbstractVector{<:BitStr}) = mean(count_vertices, samples)

"""
    soft_misloss([reg::RydbergReg], α::Real)

The soften maximum independent set loss defined as
```math
L = -1/α \\log(\\langle ψ|\\exp(α \\sum(n))|ψ\\rangle),
```
where `n` is the vertex set size.
"""
function soft_misloss(reg::RydbergReg, α::Real)
    expected = sum(zip(reg.subspace, relaxedvec(reg))) do (config, amp)
        abs2(amp) * exp(α * count_vertices(config))
    end
    -log(expected)/α
end

function logsumexp(x::AbstractArray)
    xmax = maximum(x)
    log(sum(exp.(x .- xmax))) + xmax
end

function soft_misloss(samples::AbstractVector{<:BitStr}, α::Real)
    -(logsumexp(α .* count_vertices.(samples)) - log(length(samples)))/α
end

function soft_misloss(samples::AbstractMatrix, α::Real)
    -(logsumexp(α .* sum(samples, dims=(2,))) - log(size(samples, 1)))/α
end

soft_misloss(α::Real) = reg -> soft_misloss(reg, α)
