using Plots
using RydbergEmulator
using Yao
using BenchmarkTools
using ProgressLogging
using TerminalLoggers
using Logging: global_logger

global_logger(TerminalLogger())

function benchmark_blockade(atoms::Vector{<:RydAtom})
    n = length(atoms)
    graph = unit_disk_graph(atoms, 1.0)
    s = blockade_subspace(graph)
    hs = [RydInteract(atoms, 1.0) + XTerm(n, rand()) for _ in 1:5]
    ts = rand(n)
    t = @benchmark emulate!(r, $ts, $hs) setup=(r=zero_state($n, $s))
    return minimum(t.times)
end

function benchmark_fullspace(atoms::Vector{<:RydAtom})
    n = length(atoms)
    hs = [RydInteract(atoms, 1.0) + XTerm(n, rand()) for _ in 1:5]
    ts = rand(n)
    t = @benchmark emulate!(r, $ts, $hs) setup=(r=zero_state($n))
    return minimum(t.times)
end

function run(range)
    blockade_times = Float64[]
    fullspace_times = Float64[]
    @progress for n in range
        atoms = square_lattice(n, 0.8)
        push!(blockade_times, benchmark_blockade(atoms))
        push!(fullspace_times, benchmark_fullspace(atoms))
    end
    return blockade_times, fullspace_times
end

blockade_times, fullspace_times = run(10:20)
