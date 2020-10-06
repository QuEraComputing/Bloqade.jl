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
        if n < 26
            fullspace_t = benchmark_fullspace(atoms)
            @info fullspace_t
            open("fullspace.dat", "a+") do f
                println(f, n, ", ", fullspace_t)
            end 
        end
    end
    return
end

run(10:30)

# using DelimitedFiles
# times = readdlm("benchmarks/timings.dat")
# plot(times[:, 1])
# group = repeat(["blockade", "fullspace"], inner = 11)
# xs = repeat(string.(collect(10:20)), outer = 2)

# plt = plot(legend=:topleft, yaxis="ns", xaxis="number of atoms")
# plot!(plt, string.(10:20), times[:, 1], yscale=:log10, label="blockade", markershape = :circle)
# plot!(plt, string.(10:20), times[:, 2], yscale=:log10, label="fullspace", markershape = :circle)

# groupedbar(xs, times;
#     yaxis="ns", xaxis="number of atoms",
#     yscale=:log10, bar_position=:dodge, group,
#     legend=:topleft,
# )
# savefig("benchmark.png")
# times[:,2]./times[:, 1]
