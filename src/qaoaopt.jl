export qaoa_on_graph, mean_independent_set, cmaes_train_mis

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
    mean_independent_set(graph, ϕs::AbstractVector, ts::AbstractVector; nshots=nothing)

Estimate the mean independent set by executing qaoa circuit on a graph model.
If `nshots` is nothing, there is no sampling error.
"""
function mean_independent_set(graph, ϕs::AbstractVector, ts::AbstractVector; nshots=nothing)
    reg = qaoa_on_graph(graph, ϕs, ts)
    # compute loss
    if nshots isa Nothing
        expect_mis(reg)
    else
        isets = measure_mis(reg; nshots=nshots)
        Statistics.mean(isets)
    end
end

"""
    cmaes_train_mis(graph, ϕs0, ts0)

Obtain the MIS using CMA-ES training. Return the optimal `ϕs` and `ts`.
"""
function cmaes_train_mis(graph, ϕs0, ts0)
    @assert length(ϕs0) == length(ts0)
    p = length(ϕs0)
    params = vcat(ϕs0, ts0)
    optparams = cmaes(params; num_offsprings=50, num_parents=10, maxiter=100, tol=1e-5) do params
        p = length(params)÷2
        ϕs = params[1:p]
        ts = params[p+1:end]
        -mean_independent_set(graph, ϕs, ts; nshots=nothing)
    end
    optparams[1:p], optparams[p+1:end]
end
