using Graphs
using EaRyd
using EaRyd.EaRydLattices
using Compose
using Random
using Yao

# ## The maximum independent set problem
# In graph theory, an independent set is a set of vertices in a graph, no two of which are adjacent.
# The problem of finding maximum independent sets (MIS) is NP-hard, i.e. unlikely to be solved in polynomial time.
# Even aproximating the MIS size ``\alpha(G)`` for a graph ``G=(V,E)`` is also proven to be hard.
# In this tutorial we study the MIS problem defined on diagonal-coupled unit-disk grid graphs (DUGG).
# Although these graphs have highly constraint topology, finding its MISs is still NP-hard.
# We show how to map the MIS problem on this graph to a Rydberg atom array hamiltonian,
# and use two quantum algorithms, the standard QAOA and a variational quantum algorithm with specially parametrized waveform, to find maximum independent sets.
# For those who wants to know more details, we highly recommend they to connect this tutorial with the experiment [arxiv:2202.09372](https://arxiv.org/abs/2202.09372).

# To begin with, let us create a DUGG at 0.8 filling.
# The lattice constant is 4.5μm, and blockade radius is 7.5μm.
Random.seed!(2)
atoms = generate_sites(SquareLattice(), 4, 4; scale=4.5) |> random_dropout(0.2)
img_atoms(atoms, blockade_radius=7.5)

# For the pedagogical purpose, we show the MIS size here using the graph utilities in EaRyd so that a user can compare this exact result with the quantum one.
graph = unit_disk_graph(RydAtom.(atoms), 7.5)
EaRyd.EaRydCore.exact_solve_mis(graph)

# ## QAOA with piecewise constant pulses
# The QAOA algorithm (https://arxiv.org/abs/1411.4028) can be .
# The standard definition involves applying the problem Hamiltonian ``C`` and the transverse field Hamiltonian ``B`` alternately.
# Let ``G=(V,E)`` be a graph, the hamiltonian for an MIS problem definite on it should be
# ```math
# C(G, \sigma^z) = -\sum_{i\in V} w_i \sigma_i^z + \infty \sum_{\langle i,j\rangle \in E}\sigma_i^z \sigma_j^z
# \end{equation}
# ```
# where the first summation is proportional to the size of the independent set, while the second term enfores the independence constraints.
# In a Rydberg hamiltonian, the first term is the detuning ``\Delta``.
# The second term contains an ``\infty``, which corresponds to the Rydberg blockade term that its strength decreases very fast as distance: ``\propto |r_i - r_j|^{-6}``.
# It is not a perfect independent constraint term, hence proprocessing might be required in a Rydberg atom array experiment.
#
# The transverse field Hamiltonian corresponds to the Rabi term in a Rydberg atom array.
# ```math
# B = \sum_{j=1}^{n}\sigma_j^x
# ```

# The QAOA algorithm is a classical-quantum hybrid algorithm, the classical part is the optimizer, which can be either a gradient based or non-gradient based.
# The quantum part is a Rydberg atom system evolving under a parameterized pulse sequence and finally get measured on the computational basis.
# For the convenience of simulation, we use the [`expect`](@ref) function to get the averaged measurement output and use automatic differentiation package ForwardDiff to differentiate the pulses.
# In a experimental setup, the [`expect`] should be replace by measuring on the computational basis and get the energy of bit strings averaged as the loss.
# Then one can either use non-gradient based optimizers to do the optimization or use finite difference obtain gradients of parameters.
function loss_piecewise_constant(atoms::AtomList, x::AbstractVector{T}) where T
    @assert length(x) % 2 == 0
    p = length(x)÷2
    ## detuning and rabi terms
    ## Δs = repeat(T[0.0, 1], p)
    ## Ωs = repeat(T[1.0, 0], p)
    durations = abs.(x)
    clocks = [0, cumsum(durations)...]
    Ωs = piecewise_constant(; clocks=clocks, values=repeat(T[0, 1], p))
    Δs = piecewise_constant(; clocks=clocks, values=repeat(T[1, 0], p))

    ## NOTE: check Δ
    ## hamiltonians = [rydberg_h(atoms; C=C, Ω=Ω, ϕ=zero(T), Δ=Δ) for (Δ, Ω) in zip(Δs, Ωs)]
    hamiltonian = rydberg_h(atoms; Ω=Ωs, Δ=Δs)

    subspace = blockade_subspace(atoms, 5.2)
    ## We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
    ## prob = KrylovEvolution(zero_state(Complex{T}, subspace), durations, hamiltonians)
    prob = SchrodingerProblem(zero_state(Complex{T}, subspace), sum(durations), hamiltonian)
    emulate!(prob)

    ## results are bit strings
    nbits = length(atoms)
    ## return real(sum(prob.reg.state))
    return -real(expect(sum([put(nbits, i=>ConstGate.P1) for i=1:nbits]), prob.reg))
end

# Let us check the loss function
loss_piecewise_constant(atoms, rand(4))

# The ouput shows the negative mean independent set size, where we flip the sign because in most optimizers, the loss function should be a real number to minimize.
# Here, we set the goal to maximize the average independent set size by measuring the output configuration for simplicity:
# ```math 
# \mathcal{L} = \frac{1}{N}\sum\limits_{s=1}^N n_s,
# ```
# where ``N`` is the number of measurements.
# Ideally, the goal should be maximizing the probability to get an MIS, but this value is in general not available in a real applications.
# There are some alternatives like the CVaR loss
# ```math
# \mathcal{L}_{\rm CVaR} = \frac{1}{|S|}\sum\limits_{s\sim S}n_s``, where ``S`` is the collection of best ``\alpha``(=$cvar_alpha)``N`` configurations"
# ```
# and Gibbs loss
# ```math
# \mathcal{L}_{\rm Gibbs} = -\frac{1}{\beta}\log(\langle \psi|e^{-\beta H}|\psi\rangle)
# ```
# where ``\beta`` is the "inverse temperature" as a hyperparameter.

# Let us use the forward mode automatic differentiation to get parameter gradient for free,
# for gradient based optimization.
using ForwardDiff

t0 = rand(4)

ForwardDiff.gradient(x->loss_piecewise_constant(atoms, x), t0)

# NOTE: this gradient is not consistent with the finite difference!
# It is uniform
## using FiniteDifferences

## FiniteDifferences.jacobian(central_fdm(5,1; factor=1e2), x->loss_piecewise_constant(atoms, x), t0)

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

