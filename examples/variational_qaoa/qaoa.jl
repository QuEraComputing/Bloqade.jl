# note, here we normalize parameters
function qaoa_loss(atoms::AtomList, x::AbstractVector; nshots::Int)
    @assert (length(x)+2) % 4 == 0
    npulses = (length(x)+2) ÷ 4
    # We first prepare the adiabatic pulse sequence as two piecewise linear functions
    # define the rabi waveform
    Ω_max = 2.3 * 2 * pi
    # x[1:npulses-1] is Δt for 
    clock1 = [0.0, cumsum(abs.(x[1:npulses-1]))...]
    Ω = piecewise_linear(clocks=clock1, values=x[npulses:2*npulses-1] .* Ω_max)

    # define the detuning waveform
    U = Ω_max / 2.3
    clock2 = [0.0,cumsum(abs.(x[2*npulses:3*npulses-2]))...]
    Δ = piecewise_linear(clocks=clock2, values=x[3*npulses-1:end] .* U)

    # We construct the Rydberg Hamiltonian from the defined rabi and detuning waveforms
    h = rydberg_h(atoms; C=2 * pi * 858386, Δ, Ω)

    # We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
    prob = ODEEvolution(zero_state(length(atoms)), min(clock1[end], clock2[end]), h)
    emulate!(prob)
    results = measure(prob.reg; nshots=nshots)
    -mean_rydberg(results)
end

using StochasticOptimizers
x0 = rand(14)
qaoa_loss(atoms, x0; nshots=100)

cmaes_optimizer = cmaes(x->qaoa_loss(atoms, x; nshots=100), x0; npopulation=5, noffsprings=20, σ0=0.3);
for (k, it) in enumerate(cmaes_optimizer)
    k > 100 && break
    @info "step $k, loss = $(minimum(it))"
end

minimizer(cmaes_optimizer.state)

# TODO: 
# * set correct parameters
# * more visualization
# * subspace