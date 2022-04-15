# # Solving the maximum independent set with variational optimization algorithm
# This example shows how to find the maximum independent set of a graph using adiabatic algorithm and the 
using Graphs
using Bloqade
using BloqadePlots
using Compose
using Random
using Optim         ## the optimizer

# ## The maximum independent set problem
# In graph theory, an independent set is a set of vertices in a graph, no two of which are adjacent.
# The problem of finding maximum independent sets (MIS) is NP-hard, i.e. unlikely to be solved in polynomial time.
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

# For the pedagogical purpose, we show the MIS size here using the graph utilities in Bloqade so that a user can compare this exact result with the quantum one.
# The exact MIS size can be solved using the exact solver [`exact_solve_mis`](@ref), which implements a Branching algorithm.
# It takes a graph instance as input, which can be generated by the [`unit_disk_graph`](@ref) function.
graph = BloqadeMIS.unit_disk_graph(atoms, 7.5)
BloqadeMIS.exact_solve_mis(graph)
# For advanced MIS solvers, please check [GenericTensorNetworks](@ref).

# ## The adiabatic approach
# We first prepare the adiabatic pulse sequence as two piecewise linear functions
# define the waveform for Rabi pulse
T_max = 1.65
Ω_max = 4 * 2π
Ω = piecewise_linear(clocks=[0.0, 0.2, 1.45, T_max], values=[0.0, Ω_max , Ω_max , 0])
BloqadePlots.draw(Ω)

# define the detuning waveform
Δ_start = -13 * 2π
Δ_end = 11 * 2π
Δ = piecewise_linear(clocks=[0.0, 0.2, 1.45, T_max], values=[Δ_start, Δ_start, Δ_end, Δ_end])
BloqadePlots.draw(Δ)

# Here, the total time is fixed to `T_max`, the adiabatic evolution path is specified by the piecewise linear function.
# Rydberg blockade radius can be computed with 
# ```math
# C / R_b^6 \sim \sqrt{\Delta^2 + \Omega^2}
# ```
# When `\Omega = 0`, the blockade radius for ``7.5\mu m`` is ``\sim 11 \times 2\pi MHz`` [Sepehr2022].

hamiltonian = rydberg_h(atoms; Ω=Ω, Δ=Δ)

prob = SchrodingerProblem(zero_state(nqubits(hamiltonian)), T_max, hamiltonian)
emulate!(prob)

bitstring_hist(prob.reg; nlargest=20)

using BitBasis
best_bit_string = most_probable(prob.reg, 2)

img_atoms(atoms, colors=[iszero(b) ? "white" : "black" for b in best_bit_string[1]])
img_atoms(atoms, colors=[iszero(b) ? "white" : "black" for b in best_bit_string[2]])

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

# The QAOA algorithm is a classical-quantum hybrid algorithm, the classical part is an optimizer, which can be either a gradient based or non-gradient based one.
# The quantum part is a Rydberg atom system evolving under a parameterized pulse sequence and finally get measured on the computational basis.
# For the convenience of simulation, we use the [`expect`](@ref) function to get the averaged measurement output and use automatic differentiation package ForwardDiff to differentiate the pulses.
# In a experimental setup, the [`expect`] should be replace by measuring on the computational basis and get the energy of bit strings averaged as the loss.
# Then one can either use non-gradient based optimizers to do the optimization or use finite difference obtain gradients of parameters.

# For p=3, we have
durations = [0.1, 0.5, 0.3, 0.3, 0.2, 0.4]
clocks = [0, cumsum(durations)...]
Ω2 = piecewise_constant(; clocks=clocks, values=repeat([Ω_max, 0.0], 3))
Δ2 = piecewise_constant(; clocks=clocks, values=repeat([0.0, Δ_end], 3))

hamiltonian2 = rydberg_h(atoms; Ω=Ω2, Δ=Δ2)

# piece-wise constants can be more accurately solved with the [`KrylovEvolution`](@ref) solver.
nbits = length(atoms)
prob2 = KrylovEvolution(zero_state(nbits), clocks, hamiltonian)
emulate!(prob2)

# we defined a loss function as the mean MIS size.
loss_qaoa(reg) = -real(expect(SumOfN(nsites=nbits), reg))

loss_qaoa(prob2.reg)

# This loss does not look good, we can throw it into an optimizer and see if a classical optimizer can help
# But first, let us wrap up the above code into a loss function.

function loss_piecewise_constant(atoms::AtomList, x::AbstractVector{T}) where T
    @assert length(x) % 2 == 0
    Ω_max = 4 * 2π
    Δ_end = 11 * 2π
    p = length(x)÷2
    ## detuning and rabi terms
    durations = abs.(x)
    clocks = [0, cumsum(durations)...]
    Ωs = piecewise_constant(; clocks=clocks, values=repeat(T[Ω_max, 0.0], p))
    Δs = piecewise_constant(; clocks=clocks, values=repeat(T[0.0, Δ_end], p))

    hamiltonian = rydberg_h(atoms; Ω=Ωs, Δ=Δs)

    ## We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
    subspace = blockade_subspace(atoms, 7.5)
    prob = KrylovEvolution(zero_state(Complex{T}, subspace), clocks, hamiltonian)
    emulate!(prob)

    ## results are bit strings
    nbits = length(atoms)
    ## return real(sum(prob.reg.state))
    return -real(expect(sum([put(nbits, i=>ConstGate.P1) for i=1:nbits]), prob.reg)), prob.reg
end

# !!!note
#     Run in subspace does not violate the independence constraints.
#     In practise, one needs to post-process the measured bit strings to a get a correct measure of loss.

x0 = (Random.seed!(2); rand(6))

# Let us check the loss function
mean_mis, reg0 = loss_piecewise_constant(atoms, x0)

mean_mis

# If we plot the distribution

bitstring_hist(reg0; nlargest=20)

# Let us use the non-gradient based optimizer NelderMead to optimize the loss
optresult = Optim.optimize(x->loss_piecewise_constant(atoms, x)[1], x0)

mean_mis_final, reg_final = loss_piecewise_constant(atoms, optresult.minimizer)

bitstring_hist(reg_final; nlargest=20)

# ## Smoothen Piecewise linear pulses

smooth(piecewise_linear(clocks=[0.0, 0.2, 1.45, T_max], values=[0.0, Ω_max , Ω_max , 0]); kernel_radius=0.1)

function loss_piecewise_linear(atoms::AtomList, x::AbstractVector{T}) where T
    @assert length(x) == 3
    Ω_max = 4 * 2π
    Δ_start = -13 * 2π
    Δ_end = 11 * 2π
    Δ0 = 11 * 2π
    T_max = 0.6
    ## detuning and rabi terms
    Δs = smooth(piecewise_linear(clocks=T[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, T_max], values=T[Δ_start, Δ_start, Δ0*x[1], Δ0*x[2], Δ0*x[3], Δ_end, Δ_end]); kernel_radius=0.1)
    Ωs = smooth(piecewise_linear(clocks=T[0.0, 0.1, 0.5, T_max], values=T[0.0, Ω_max , Ω_max , 0]); kernel_radius=0.05)

    hamiltonian = rydberg_h(atoms; Ω=Ωs, Δ=Δs)

    ## We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
    subspace = blockade_subspace(atoms, 7.5)
    prob = SchrodingerProblem(zero_state(Complex{T}, subspace), T_max, hamiltonian)
    emulate!(prob)

    ## results are bit strings
    nbits = length(atoms)
    ## return real(sum(prob.reg.state))
    return -real(expect(sum([put(nbits, i=>ConstGate.P1) for i=1:nbits]), prob.reg)), prob.reg, Δs
end

x0 = [0.1, 0.2, 0.2]

# Let us check the loss function
mean_mis, reg0, Δ_initial = loss_piecewise_linear(atoms, x0)
BloqadePlots.draw(Δ_initial)

mean_mis

# If we plot the distribution

bitstring_hist(reg0; nlargest=20)

# Let us use the non-gradient based optimizer NelderMead to optimize the loss
optresult = Optim.optimize(x->loss_piecewise_linear(atoms, x)[1], x0)

mean_mis_final, reg_final, Δ_final = loss_piecewise_linear(atoms, optresult.minimizer)

bitstring_hist(reg_final; nlargest=20)

BloqadePlots.draw(Δ_final)

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
