using Plots
using StatsPlots
using RydbergEmulator
using Yao
using BenchmarkTools
using DelimitedFiles

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
    for n in range
        atoms = square_lattice(n, 0.8)
        @info n
        blockade_t = benchmark_blockade(atoms)
        @info blockade_t
        open("blockade.dat", "a+") do f
            println(f, n, ", ", blockade_t)
        end
    end
    return
end

run(26:40)

using DelimitedFiles
using Unitful: ns, ms
blockade_times = readdlm("benchmarks/blockade.dat")
fullspace_times = readdlm("benchmarks/fullspace.dat")

to_ms(x) = map(x->x.val, uconvert.(ms, x .* ns))

plt = plot(
    legend=:topleft,
    xticks=9:2:37,
    yaxis="ms", xaxis="number of atoms",
    size=(1000, 600),
    xtickfontsize=14, ytickfontsize=14,
    xguidefontsize=14, yguidefontsize=14,
    legendfont=14,
    thickness_scaling=1.5,
)
plot!(plt, 10:36, to_ms(blockade_times[:, 2]), yscale=:log10, label="blockade", markershape = :circle)
plot!(plt, 10:24, to_ms(fullspace_times[:, 2]), yscale=:log10, label="fullspace", markershape = :circle)

savefig("benchmark.png")
