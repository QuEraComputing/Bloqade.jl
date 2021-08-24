"""
    is_independent_set(graph::AbstractGraph, config)

Return `true` if `config` is an independent set of graph.
`config` can be a `BitStr`, a vector, or any iterable.
"""
function is_independent_set(graph::AbstractGraph, config)
    for i in 1:nv(graph)
        for j in i+1:nv(graph)
            config[i] == 1 && config[j] == 1 && has_edge(graph, i, j) && return false
        end
    end
    return true
end

"""
    to_independent_set!(graph::AbstractGraph, config::AbstractVector)

Eliminate vertices in `config` so that remaining vertices do not have connected edges.
This algorithm is a naive vertex elimination that does not nesesarily give the maximum possible vertex set.

```@example
# run the following code in Atom/VSCode
atoms = RydAtom.([(0.0, 1.0), (1.0, 0.), (2.0, 0.0),
    (1.0, 1.0), (1.0, 2.0), (2.0, 2.0)])
graph = unit_disk_graph(atoms, 1.5)

config = [1, 1, 1, 0, 1, 1]
viz_config(atoms, graph, config)

to_independent_set!(graph, config)
viz_config(atoms, graph, config)
```
"""
function to_independent_set!(graph::AbstractGraph, config::AbstractVector)
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

function num_mis_violation(config, graph::AbstractGraph, i::Int)
    ci = config[i]
    ci == 1 || return 0
    count(j -> config[j] == 1 && has_edge(graph, i, j), 1:length(config))
end

"""
    exact_solve_mis(g::AbstractGraph)

Return the exact MIS size of a graph `g`.
"""
exact_solve_mis(g::AbstractGraph) = mis2(EliminateGraph(adjacency_matrix(g)))

"""
    reduced_mean_rydberg(reg::RydbergReg, graph::AbstractGraph; add_vertices::Bool=false)

independent set eliminated mean rydberg.
"""
function reduced_mean_rydberg(reg::RydbergReg, graph::AbstractGraph; add_vertices::Bool=false)
    mean_ryd = zero(real(eltype(reg.state)))
    for (c, amp) in zip(vec(reg.subspace), vec(reg.state))
        new_c = to_independent_set!(graph, bitarray(c, nv(graph)))
        add_vertices && add_vertices!(new_c, graph)
        mean_ryd += abs2(amp) * count_vertices(packbits(new_c))
    end
    return mean_ryd
end

function add_vertices!(config::AbstractVector, graph::AbstractGraph)
    for (k, c) in enumerate(config)
        if iszero(c)
            config[k] = 1
            if is_independent_set(graph, config)
                continue
            else
                config[k] = 0
            end
        end
    end
    return config
end

"""
    independent_set_probabilities(reg::RydbergReg, graph::AbstractGraph, mis::Int = exact_solve_mis(graph); add_vertices::Bool = false)

Return a `Vector` of the probabilities of all independent set size.
An optional argument `mis` can be applied to use pre-calculated mis size.

The probability of each size can be queried via `prob[size+1]`, e.g

```julia
prob = independent_set_probabilities(r, graph)
prob[end]   # the probability of mis
prob[end-1] # the probability of mis - 1
```
"""
function independent_set_probabilities(reg::RydbergReg, graph::AbstractGraph, mis::Int = exact_solve_mis(graph); add_vertices::Bool = false)
    probs = zeros(real(eltype(reg.state)), mis+1)

    for (c, amp) in zip(vec(reg.subspace), vec(reg.state))
        new_c_array = to_independent_set!(graph, bitarray(c, nv(graph)))
        add_vertices && add_vertices!(new_c_array, graph)
        new_c = packbits(new_c_array)

        nvertices = count_vertices(new_c)
        probs[nvertices+1] += abs2(amp)
    end
    return probs
end
