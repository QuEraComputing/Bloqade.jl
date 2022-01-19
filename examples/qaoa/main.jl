using Graphs, EaRyd, EaRyd.EaRydLattices
using Compose

function vizconfig(atoms::AtomList; config)
    colors=map(c->iszero(c) ? "blue" : "red", config)
    img, (dx, dy) = EaRydLattices.img_atoms(atoms; scale=1.5, colors=colors)
    img
end

# create a diagonal coupled square lattice with 0.7 filling.
atoms = generate_sites(SquareLattice(), 4, 4; scale=9.629) |> random_dropout(0.3)
vizconfig(atoms, config=rand(Bool, length(atoms)))

# We first prepare the adiabatic pulse sequence as two piecewise linear functions
# define the rabi waveform
Ω_max = 2.3 * 2 * pi
Ω = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[0.0, Ω_max , Ω_max , 0])

# define the detuning waveform
U = Ω_max / 2.3
Δ = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[-6*U, -6*U, 2*U , 2*U])

# We construct the Rydberg Hamiltonian from the defined rabi and detuning waveforms
h = rydberg_h(atoms; C=2 * pi * 858386, Δ, Ω)

# We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
prob = ODEEvolution(zero_state(length(atoms)), 1.6, h)
emulate!(prob)
results = measure(prob.reg; nshots=100)

# negative mean loss
using BitBasis
function negative_mean_mis(spinconfigs::AbstractVector{<:BitStr})
    -sum(count_ones, spinconfigs)/length(spinconfigs)
end

negative_mean_mis(results)

# note, here we normalize parameters
function qaoa_loss(atoms::AtomList, x::AbstractVector; nshots::Int)
    @assert (length(x)+2) % 4 == 0
    npulses = (length(x)+2) ÷ 4
    # We first prepare the adiabatic pulse sequence as two piecewise linear functions
    # define the rabi waveform
    Ω_max = 2.3 * 2 * pi
    clock1 = [0.0,cumsum(abs.(x[1:npulses-1]))...]
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
    negative_mean_mis(results)
end

using StochasticOptimizers
x0 = rand(14)

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
