using CSV
using DataFrames

using Plots
using LaTeXStrings

using Bloqade
using BloqadeNoisy

source = Base.source_dir()
data_dir = source*"/data/"
figure_dir = source*"/figures/"

#Example 1: short 15MHz resonant Rabi

whitepaper_data = CSV.read(data_dir*"15MHz.csv", DataFrame, delim = ",", header = false)

times = collect(whitepaper_data[1, :])
data = collect(whitepaper_data[2, :])

reg = zero_state(1)
h = rydberg_h([(0,0)], Ω = 15.7, Δ = 0) #Ω=15 has some overall bias which the simulation does not account for
save_times = LinRange(0, last(times), 500)

ns = NoisySchrodingerProblem(reg, save_times, h, Aquila())
sim = emulate(ns, 1000, [mat(Op.n)]; readout_error = true, report_error = true)

scatter(times, data, label = "whitepaper", xlabel = "t (µs)", ylabel = L"\langle n \rangle", color = :red, legendfontsize = 10, tickfontsize = 12, labelfontsize = 14, legend = :outerbottomright)
plot!(save_times, sim.expectations[1], ribbon = sim.propagated_error[1], label = "simulated", color = :blue)
xlims!(first(times), last(times))
png(figure_dir*"15MHz_short.png")

#Example 2: long 15MHz resonant Rabi

whitepaper_data = CSV.read(data_dir*"15MHz_long.csv", DataFrame, delim = ",", header = false)

times = collect(whitepaper_data[1, :])
data = collect(whitepaper_data[2, :])

reg = zero_state(1)
h = rydberg_h([(0,0)], Ω = 15.3, Δ = 0) #Ω is not known exactly again
save_times = LinRange(0, last(times), 500)

ns = NoisySchrodingerProblem(reg, save_times, h, Aquila())
sim = emulate(ns, 2000, [mat(Op.n)]; readout_error = true)

plot(times, data, label = "whitepaper", xlabel = "t (µs)", ylabel = L"\langle n \rangle", color = :red, legendfontsize = 10, tickfontsize = 12, labelfontsize = 14)
plot!(save_times, sim, label = "simulated", color = :blue)
png(figure_dir*"15MHz_long_incoherent.png")

#Example 3: long detuned Rabi

whitepaper_data = CSV.read(data_dir*"15MHz_detuned_long.csv", DataFrame, delim = ",", header = false)

times = collect(whitepaper_data[1, :])
data = collect(whitepaper_data[2, :])

reg = zero_state(1)
h = rydberg_h([(0,0)], Ω = 14.5, Δ = 17.5)
save_times = LinRange(0, last(times), 500)
ns = NoisySchrodingerProblem(reg, save_times, h, Aquila())

sim = emulate(ns, 2000, [mat(Op.n)]; readout_error = true)

plot(times, data, label = "whitepaper", xlabel = "t (µs)", ylabel = L"\langle n \rangle", color = :red, legendfontsize = 10, tickfontsize = 12, labelfontsize = 14)
plot!(save_times, sim, label = "simulated", color = :blue)
png(figure_dir*"15MHz_detuned_long.png")

#Example 4: 3-atom Blockade-enhanced Rabi

whitepaper_data = CSV.read(data_dir*"3q_blockaded_rabi.csv", DataFrame, delim = ",", header = false)

times = collect(whitepaper_data[1, :])
data = collect(whitepaper_data[2, :])

atoms = [ #From whitepaper
            (
              0.0,
              4e-06
            ),
            (
              -2.6e-06,
              0.0
            ),
            (
              2.6e-06,
              0.0
            ),
]
atoms = [a .* 1f6 for a in atoms] #units

reg = zero_state(3)
h = rydberg_h(atoms, Ω = 2.58*2π, Δ = 0)
save_times = LinRange(0, last(times), 500)

ns = NoisySchrodingerProblem(reg, save_times, h, Aquila())

cmat = Aquila().confusion_mat(3) #get confusion mat to calculate ground state probability
sim = emulate(ns, 1000, sol -> [(cmat*abs.(u).^2)[1] for u in sol])

plot(times, data, label = "whitepaper", xlabel = "t (µs)", ylabel = L"p_{000}",color = :red, legendfontsize = 10, tickfontsize = 12, labelfontsize = 14)
plot!(save_times, simulation_series_mean(sim), label = "simulated", color = :blue)
png(figure_dir*"3_atom_blockaded.png")