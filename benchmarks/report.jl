using Plots
using RydbergEmulator
using Yao
using BenchmarkTools

function benchmark_blockade(atoms)
    graph = unit_disk_graph(atoms, 1.0)
    s = blockade_subspace(graph)
    hs = [RydInteract(atoms, 1.0) + XTerm(n, rand()) for _ in 1:5]
    ts = rand(n)
    t = @benchmark emulate!(r, $ts, hs) setup=(r=zero_state($n, $s))
    return minimum(t.times)
end

function benchmark_fullspace(atoms)
    hs = [RydInteract(atoms, 1.0) + XTerm(n, rand()) for _ in 1:5]
    ts = rand(n)
    t = @benchmark emulate!(r, $ts, hs) setup=(r=zero_state($n))
    return minimum(t.times)
end

function run(range)
    times = Float64[]
    for n in range
        atoms = square_lattice(n, 0.8)
        t = benchmark_blockade(atoms)
        push!(times, t)
    end
    return times
end

times = run(10:20)
