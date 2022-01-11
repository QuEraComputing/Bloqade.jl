using EaRyd
using Random
Random.seed!(42)
# build lattice
nsites = 10
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)
h = rydberg_h(atoms;C = 2π * 858386, Ω=4π)
config = rand(0:1, 10)
init = product_state(config)

# # construct time evolution
iteration = 1:120
ts = [0.01 for _ in iteration];
hs = [h for _ in iteration];
prob = KrylovEvolution(init, ts, hs)

# measuring obervable
clocks = cumsum(ts)
entropy = zeros(length(iteration))
domain_mat = zeros(nsites-1, length(iteration))
density_mat = zeros(nsites, length(iteration))

for info in prob
    for i in 1:nsites
        density_mat[i, info.step] = expect(put(nsites, i=>Op.n), info.reg)
    end

    for i in 1:nsites-1
        corr = real(expect(put(nsites, (i, i+1)=>kron(Op.n, Op.n)), info.reg))
        obs = density_mat[i, info.step] + density_mat[i+1, info.step] - 2corr
        domain_mat[i, info.step] = obs
    end

    rho = density_matrix(info.reg, (1,2,3,4,5))
    entropy[info.step] = von_neumann_entropy(rho)
end

using CairoMakie
fig = Figure(size=(10, 5));
ax = Axis(fig[1, 1])

for i in 1:nsites
    lines!(clocks, density_mat[i, :])
end
fig

heatmap(clocks, 1:nsites, density_mat')
heatmap(clocks, 1:nsites, domain_mat')
domain_avg = vec(sum(domain_mat, dims=1)/(nsites-1))
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, clocks, domain_avg)
lines!(ax, clocks, entropy)
fig


# TODO: subspace

