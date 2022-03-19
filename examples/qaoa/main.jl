using Graphs
using EaRyd
using EaRyd.EaRydLattices
using Compose
using Random
using Yao

# create a diagonal coupled square lattice with 0.8 filling.
# In the experiment [arxiv:2202.09372](https://arxiv.org/abs/2202.09372),
# The lattice constant is 4.5μm, and blockade radius is 7.5μm.
Random.seed!(2)
atoms = generate_sites(SquareLattice(), 4, 4; scale=4.5) |> random_dropout(0.2)

# ## Piecewise constant pulses
function loss_piecewise_constant(atoms::AtomList, x::AbstractVector{T}) where T
    @assert length(x) % 2 == 0
    ## durations can not be negative
    durations = ones(length(x))

    ## detuning and rabi terms
    #Δs = repeat(T[0.0, 1], length(x)÷2)
    #Ωs = repeat(T[1.0, 0], length(x)÷2)
    Ωs = piecewise_constant(; clocks=[0, cumsum(durations)...], values=vcat([[x[2i-1], zero(T)] for i in 1:length(x)÷2]...))
    Δs = piecewise_constant(; clocks=[0, cumsum(durations)...], values=vcat([[zero(T), x[2i]] for i in 1:length(x)÷2]...))

    ## hamiltonians at different time step
    C = T(2π * 858386)

    ## NOTE: check Δ
    hamiltonians = [rydberg_h(atoms; C=C, Ω=Ω, ϕ=zero(T), Δ=Δ) for (Δ, Ω) in zip(Δs, Ωs)]
    #hamiltonian = rydberg_h(atoms; C=C, Ω=Ωs, Δ=Δs)

    subspace = blockade_subspace(atoms, 5.2)
    ## We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
    prob = KrylovEvolution(zero_state(Complex{T}, subspace), durations, hamiltonians)
    #prob = ODEEvolution(zero_state(Complex{T}, subspace), sum(durations), hamiltonian)
    emulate!(prob)

    ## results are bit strings
    nbits = length(atoms)
    return real(expect(sum([put(nbits, i=>ConstGate.P1) for i=1:nbits]), prob.reg))
end

using ForwardDiff

ForwardDiff.gradient(x->loss_piecewise_constant(atoms, x), rand(4))

# Let us use the forward mode automatic differentiation to get parameter gradient for free,
# for gradient based optimization.

# We first prepare the adiabatic pulse sequence as two piecewise linear functions
# define the rabi waveform
Ω_max = 2.3 * 2 * pi
Ω = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[0.0, Ω_max , Ω_max , 0])

# define the detuning waveform
U = Ω_max / 2.3
Δ = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[-6*U, -6*U, 2*U , 2*U])

# We construct the Rydberg Hamiltonian from the defined rabi and detuning waveforms
# check C?
# smoothen the curve?
h = rydberg_h(atoms; Δ, Ω)

# We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
prob = ODEEvolution(zero_state(length(atoms)), 1.6, h)
emulate!(prob)
results = measure(prob.reg; nshots=100)

# add a new API for general Vector
graph = unit_disk_graph(RydAtom.(atoms), 7.5)
EaRyd.EaRydCore.exact_solve_mis(graph)

function sample_piecewise_constant(atoms::AtomList, x::AbstractVector; nshots::Int, blockade_radius=5.2)
    @assert length(x) % 3 == 0
    npulses = length(x) ÷ 3
    # We first prepare the adiabatic pulse sequence as two piecewise linear functions
    # define the rabi waveform
    Ω = 4.0

    durations = abs.(x[1:npulses])

    # phases
    ϕs = x[2*npulses+1:3*npulses] .* 2π

    # hamiltonians at different time step
    C = 2π * 858386

    # the detuning
    Δ = 4π

    # NOTE: check Δ
    hamiltonians = [rydberg_h(atoms; C=C, Ω=Ω, ϕ=ϕ, Δ=Δ) for ϕ in ϕs]

    subspace = blockade_subspace(atoms, blockade_radius)
    # We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
    prob = KrylovEvolution(zero_state(subspace), durations, hamiltonians)
    emulate!(prob)

    # results are bit strings
    return measure(prob.reg; nshots=nshots)
end

function qaoa_loss_piecewise_constant(atoms::AtomList, x::AbstractVector; nshots::Int)
    bitstrings = sample_piecewise_constant(atoms, x; nshots=nshots)
    # call the loss function
    return -mean_rydberg(bitstrings)
end

x0 = rand(9)
qaoa_loss_piecewise_constant(atoms, x0; nshots=100)

using Optim
optimize_result = Optim.optimize(x->qaoa_loss_piecewise_constant(atoms, x; nshots=100), x0, NelderMead(), Optim.Options(show_trace=true, iterations=100))
println("The final loss is $(minimum(optimize_result))")

# sample some configurations
configs = sample_piecewise_constant(atoms, Optim.minimizer(optimize_result); nshots=100)
masks = [EaRydCore.is_independent_set(config, graph) for config in configs]
valid_configs = configs[masks]
missize, best_index = findmax(count_ones, valid_configs)
best_config = valid_configs[best_index]
println("The best MIS is $(best_config), its MIS size is $missize")
vizconfig(atoms, config=collect(best_config))

