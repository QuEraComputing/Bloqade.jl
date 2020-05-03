export qaoa_on_graph, mean_independent_set

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
        sum(isets)/length(isets)
    end
end
