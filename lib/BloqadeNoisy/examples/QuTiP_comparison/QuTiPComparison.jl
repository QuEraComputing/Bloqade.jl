using Bloqade
using BloqadeNoisy

using LaTeXStrings
using Plots
using CSV
using DataFrames

reg = zero_state(3)
h = rydberg_h([(0.0, 0.0), (8.0, 0.0), (18.0, 0.0)], Ω=2π, Δ=0)
tend = 4.0
save_times = LinRange(0, tend, 200)

#Add collapse operators to simulation
c_ops = [mat(put(3, i => sqrt(0.1) * (X + Y * im) / 2)) for i in 1:3]
ns = NoisySchrodingerProblem(reg, save_times, h, c_ops)

#get noisy expectation values
ntraj = 2000
expecs = [mat(put(3, i => Op.n)) for i in 1:3]
@time sim = emulate(ns, ntraj, expecs; report_error = true, ensemble_algo=EnsembleSerial())

plot(xlabel="t (µs)", ylabel=L"\langle n_i \rangle")
for (i, _) in enumerate(expecs)
    plot!(save_times,
        sim.expectations[i],
        ribbon=sim.twosigma[i],
        color=:blue,
        label=(i == 1 ? "BloqadeNoisy" : "")
    )
end

qutip_vals = CSV.read("/Users/queraintern/Documents/GitHub/Bloqade.jl/lib/BloqadeNoisy/examples/QuTiP_comparison/3q_expect.csv", DataFrame, header=false, delim=",")
qutip_vals_noiseless = CSV.read("/Users/queraintern/Documents/GitHub/Bloqade.jl/lib/BloqadeNoisy/examples/QuTiP_comparison/3q_expect_nonoise.csv", DataFrame, header=false, delim=",")

for i in 1:3
    plot!(
        0:1.0f-2:3.99,
        collect(qutip_vals[i, :]),
        color=:red,
        linestyle=:dash,
        label=i == 1 ? "QuTiP" : ""
    )
end

for i in 1:3
    plot!(
        0:1.0f-2:3.99,
        collect(qutip_vals_noiseless[i, :]),
        color=:black,
        label=i == 1 ? "noiseless" : "",
        alpha=0.3
    )
end

current()
ylims!(0, 1.1)

png("qutip_comparison.png")