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
using [`add_random_vertices`](@ref).

```julia
mean_rydberg(r) do config
    config = to_independent_set(config, graph)
    add_random_vertices(config, graph, 10)
    return config
end
```

Or one can also just add vertice by atom order

```julia
mean_rydberg(r) do config
    config = to_independent_set(config, graph)
    add_vertices!(config, graph)
    return config
end
```
"""
function mean_rydberg end

mean_rydberg(x) = mean_rydberg(identity, x)

struct SubspaceMap
    d::Dict{Int, Int}
end

function SubspaceMap(f, subspace::Subspace)
    d = Dict{Int, Int}()
    sizehint!(d, length(subspace))
    origin = vec(subspace)
    Threads.@threads for cfg in origin
        d[cfg] = to_int64(f(cfg))
    end
    return SubspaceMap(d)
end

Base.length(map::SubspaceMap) = length(map.d)
Base.getindex(map::SubspaceMap, cfg::Int) = map.d[cfg]
(map::SubspaceMap)(cfg::Int) = map[cfg]

function to_int64(b::BitVector) # workaround type piracy
    length(b) <= 64 || throw(ArgumentError("length is larger than 64"))
    # NOTE: since we only use this to calculate number of ones
    # thus we don't need to check top bit
    return reinterpret(Int, b.chunks[1])
end

struct ConfigAmplitude{Reg <: Yao.AbstractRegister}
    reg::Reg
    range::UnitRange{Int}
end

ConfigAmplitude(reg::Yao.AbstractRegister{1}) = ConfigAmplitude(reg, 1:size(reg.state, 1))

Base.eltype(it::ConfigAmplitude) = Tuple{Int, Yao.datatype(it.reg)}
Base.length(it::ConfigAmplitude) = length(it.range)

function Transducers.halve(it::ConfigAmplitude)
    left, right = Transducers.halve(it.range)
    return ConfigAmplitude(it.reg, left), ConfigAmplitude(it.reg, right)
end

function Base.iterate(it::ConfigAmplitude{<:RydbergReg{ComplexLayout}}, idx::Int = first(it.range))
    idx > last(it.range) && return
    cfg = it.reg.subspace.subspace_v[idx]
    @inbounds amp = it.reg.state[idx]
    return (cfg, amp), idx + 1
end

function Base.iterate(it::ConfigAmplitude{<:RydbergReg{RealLayout}}, idx::Int = first(it.range))
    idx > last(it.range) && return
    cfg = it.reg.subspace.subspace_v[idx]
    @inbounds amp = it.reg.state[idx, 1] + im * it.reg.state[idx, 2]
    return (cfg, amp), idx + 1
end

function Base.iterate(it::ConfigAmplitude{<:Yao.ArrayReg}, idx::Int = first(it.range))
    idx > last(it.range) && return
    return (idx-1, it.reg.state[idx]), idx + 1
end

function mean_rydberg(f, reg::Yao.AbstractRegister)
    return ThreadsX.sum(ConfigAmplitude(reg)) do (c, amp)
        nvertices = count_vertices(f(c))
        return abs2(amp) * nvertices
    end
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
    expected = ThreadsX.sum(ConfigAmplitude(reg)) do (config, amp)
        abs2(amp) * exp(α * count_vertices(f(config)))
    end
    return -log(expected) / α 
end

function logsumexp(x::AbstractArray)
    xmax = maximum(x)
    return log(sum(exp.(x .- xmax))) + xmax
end

function gibbs_loss(f, samples::AbstractVector, α::Real)
    expect = ThreadsX.map(x->α * count_vertices(f(x)), samples)
    return -(logsumexp(expect) - log(length(samples)))/α
end

function gibbs_loss(f, samples::AbstractMatrix, α::Real)
    -(logsumexp(α .* sum(samples, dims=(2,))) - log(size(samples, 1)))/α
end

gibbs_loss(α::Real) = reg -> gibbs_loss(reg, α)

"""
    is_independent_set(config, graph::AbstractGraph)

Return `true` if `config` is an independent set of graph.
`config` can be a `BitStr`, a vector, or any iterable.
"""
function is_independent_set(config, graph::AbstractGraph)
    for i in 1:nv(graph)
        for j in i+1:nv(graph)
            config[i] == 1 && config[j] == 1 && has_edge(graph, i, j) && return false
        end
    end
    return true
end

"""
    to_independent_set!(config::AbstractVector, graph::AbstractGraph)

Eliminate vertices in `config` so that remaining vertices do not have connected edges.
This algorithm is a naive vertex elimination that does not nesesarily give the maximum possible vertex set.

```@example
# run the following code in Atom/VSCode
atoms = RydAtom.([(0.0, 1.0), (1.0, 0.), (2.0, 0.0),
    (1.0, 1.0), (1.0, 2.0), (2.0, 2.0)])
graph = unit_disk_graph(atoms, 1.5)

config = [1, 1, 1, 0, 1, 1]
viz_config(atoms, graph, config)

to_independent_set!(config, graph)
viz_config(atoms, graph, config)
```
"""
function to_independent_set!(config::AbstractVector, graph::AbstractGraph)
    N = length(config)
    n = typemax(Int)
    while true
        n, loc = findmax(map(i->num_mis_violation(config, graph, i), 1:N))
        if n == 0
            break
        else
            config[loc] = 0
        end
    end
    return config
end

function to_independent_set(config::Integer, graph::AbstractGraph)
    return to_independent_set!(bitarray(config, nv(graph)), graph)
end

"""
    num_mis_violation(config, graph::AbstractGraph, i::Int)

Calculate the number of MIS violations for `i`-th vertex in `graph`
and configuration `config`. The `config` should be a subtype of
`AbstractVector`.
"""
Base.@propagate_inbounds function num_mis_violation(config, graph::AbstractGraph, i::Int)
    ci = config[i]
    ci == 1 || return 0

    return count(1:length(config)) do j
        config[j] == 1 && has_edge(graph, i, j)
    end
end

"""
    exact_solve_mis(g::AbstractGraph)

Return the exact MIS size of a graph `g`.
"""
exact_solve_mis(g::AbstractGraph) = mis2(EliminateGraph(adjacency_matrix(g)))

"""
    add_random_vertices([rng=GLOBAL_RNG], config::AbstractVector, graph::AbstractGraph, ntrials::Int = 10)

Add vertices randomly to given configuration for `ntrials` times and
pick the one that has largest [`count_vertices`](@ref).

# Arguments

- `rng`: optional, Random Number Generator.
- `config`: configuration to tweak.
- `graph`: problem graph.
- `ntrials`: number of trials to use, default is `10`.
"""
function add_random_vertices(
    config::AbstractVector,
    graph::AbstractGraph,
    ntrials::Int = 10)

    return add_random_vertices(
        Random.GLOBAL_RNG,
        config,
        graph,
        ntrials
    )
end

function add_random_vertices(
        rng::AbstractRNG,
        config::AbstractVector,
        graph::AbstractGraph,
        ntrials::Int = 10
    )
    perm = collect(1:nv(graph))
    nvertices = count_vertices(config)
    origin = config
    maximum_config = copy(origin)
    config_candidate = copy(origin)

    for _ in 1:ntrials
        randperm!(rng, perm)
        copyto!(config_candidate, origin)
        add_vertices!(config_candidate, graph, perm)
        nvertices_candidate = count_vertices(config_candidate)
        if nvertices_candidate > nvertices
            copyto!(maximum_config, config_candidate)
            nvertices = nvertices_candidate
        end
    end
    return maximum_config
end

add_vertices(config::AbstractVector, graph::AbstractGraph, perm=eachindex(config)) =
    add_vertices!(copy(config), graph, perm)

function add_vertices!(config::AbstractVector, graph::AbstractGraph, perm=eachindex(config))
    for (k, c) in zip(perm, config)
        if iszero(c)
            config[k] = 1
            if is_independent_set(config, graph)
                continue
            else
                config[k] = 0
            end
        end
    end
    return config
end

"""
    independent_set_probabilities([f], reg::Yao.AbstractRegister, graph_or_mis)

Calculate the probabilities of independent sets with given postprocessing function
`f(config) -> config`. The default postprocessing function `f` will only reduce all
configurations to independent set.

# Arguments

- `f`: optional, postprocessing function, default is [`to_independent_set`](@ref).
- `reg`: required, the register object.
- `graph_or_mis`: a problem graph or the MIS size of the problem
    graph (can be calculated via [`exact_solve_mis`](@ref)).
"""
function independent_set_probabilities end

function independent_set_probabilities(reg::Yao.AbstractRegister, graph::AbstractGraph)
    independent_set_probabilities(reg, graph) do config
        to_independent_set(config, graph)
    end
end

function independent_set_probabilities(f, reg::Yao.AbstractRegister, graph::AbstractGraph)
    return independent_set_probabilities(f, reg, exact_solve_mis(graph))
end

function independent_set_probabilities(f, reg::Yao.AbstractRegister, mis::Int)
    v2amp = ThreadsX.map(ConfigAmplitude(reg)) do (c, amp)
        return count_vertices(f(c)), amp
    end

    return ThreadsX.map(0:mis) do k
        ThreadsX.sum(v2amp) do (nvertices, amp)
            if nvertices == k
                abs2(amp)
            else
                zero(real(typeof(amp)))
            end
        end
    end
end

"""
    mis_postprocessing(config, graph::AbstractGraph; ntrials::Int=10)

The postprocessing protocal used in Harvard.

# Arguments

- `config`: configuration to postprocess.
- `graph`: the problem graph.

# Keyword Arguments

- `ntrials`: number of trials to use.
"""
function mis_postprocessing(config, graph::AbstractGraph; ntrials::Int=10)
    config = to_independent_set(config, graph)
    config = add_random_vertices(config, graph, ntrials)
    return config
end

"""
    mis_postprocessing(graph::AbstractGraph; ntrials::Int = 10)

Curried version of [`mis_postprocessing`](@ref).

# Example

to calculate `mean_rydberg` loss with postprocessing used in
Harvard experiment.

```julia
mean_rydberg(mis_postprocessing(graph), reg)
```
"""
function mis_postprocessing(graph::AbstractGraph; ntrials::Int = 10)
    return function postprocessing(config)
        mis_postprocessing(config, graph; ntrials)
    end
end
