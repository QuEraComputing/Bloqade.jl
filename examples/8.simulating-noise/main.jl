# # [Simulating Noise](@id blockade)

# ## Background
#
# ### The Lindblad master equation
# Open quantum systems in the limit of ultraweak coupling to a Markovian bath can be modelled using the Lindblad master equation:
# ```math
# \frac{\partial \rho}{\partial t} = i[\rho, \mathcal H_{eff}] + \sum_{k}\gamma_k L_k\rho L_k^\dagger
# ``
# where `` \mathcal H_{eff}`` is the effective Hamiltonian, `` L_k`` are the quantum jump operators, and `` \gamma_k`` are the jump rates. The jump operators describe the coupling to the bath, and they are without loss of generality taken to be traceless. The effective Hamiltonian is non-Hermitian and is related to the closed-system Hamiltonian ``\mathcal H`` via ``\mathcal H_{eff} = \mathcal H-\frac{i}{2}\sum_{k}\gamma_k L_k^\dagger L_k``. 
#
# ### Stochastic wavefunction method
# The infinitesimal form of this channel can be put into Kraus map form as
# ```math
# \rho(t+dt) = (1-i dt \mathcal H_{eff})\rho(t) (1+i dt \mathcal H_{eff}) + dt\sum_{k}\gamma_k L_k \rho(t) L_k^\dagger
# ``
# This corresponds to a quantum jump `` L_k`` with probability `` dp_k = dt\gamma_k\operatorname{Tr}(L_k \rho L_k^\dagger)``. If `` \rho = |\psi\rangle\langle \psi|``is a (normalized) pure state, then 
# `` dp_k = dt \gamma_k \Vert L_k|\psi\rangle\Vert^2``, the norm of the state after undergoing the quantum jump. With probability `` 1-dp`` where `` dp = \sum_k dp_k`` is the total probability of experencing a quantum jump, the system evolves to ``(1-i dt \mathcal H_{eff})\rho (1+i dt \mathcal H_{eff}) \approx \rho + i dt[\rho, \mathcal H_{eff}]``. This corresponds to the normal Liouville-Von Neumann equation with the non-Hermitian effective Hamiltonian ``\mathcal H_{eff}``. The physical interpretation of this is that an absence of a quantum jump also has an affect on the system evolution.
#
# These dynamics can be modeled stochastically by chopping the evolution into small intervals `` dt``. At each time step, the state `` |\psi\rangle`` evolves to `` \frac{1}{\mathcal N}(1-idt\mathcal H_{eff})|\psi \rangle``
# with probability ``1-dp`` and evolves to ``\frac{1}{\mathcal N}L_k |\psi \rangle`` with probability ``dp_k``, where ``\frac{1}{\mathcal N}``  is the appropriate normalization. Performing this stochastic evolution results in a single trajectory consisting of a number of random quantum jumps at random times. Averaging over ``j \in \{1,...,m\}`` trajectories results in an ensemble of state ``|\psi_j(t)\rangle``, and averaging ``|\psi_j\rangle\langle\psi_j|`` over these trajectories converges to the density matrix ``\rho`` produced by the master equation in the limit as ``m\to \infty``.
#
# ### Stochastic wavefunction algorithm
# In practice, directly simulating the stochastic process described above is computationally inefficient. Noticing that 
# ```math
# \Vert(1-idt\mathcal H_{eff})|\psi \rangle\Vert^2 = 1 - dt\sum_k \gamma_k\langle \psi| L_k^\dagger L_k |\psi\rangle = 1 - dp \ ,
# ```
# the change in the wavefunction norm over an interval of coherent evolution is related to the accumulation of probability `` p`` of experiencing a quantum jump. The result is the following algorithm:
# 1. Choose a random number `` r``
# 2. Evolve according to `` \frac{\partial}{\partial t}|\psi\rangle = -i\mathcal H_{eff}|\psi\rangle`` without normalizing the state
# 3. When `` \Vert |\psi\rangle \Vert^2 = r``, trigger a quantum jump `` k`` according to the distribution `` p_k \equiv dp_k/\sum_k dp_k`` and normalize
#
# ### Further reading
# For a more in-depth description of the algorithm, please refer to refs. [1,2].

using BloqadeNoisy

using Bloqade
using Yao
using CSV
using DataFrames
using JSON
using LaTeXStrings
using LinearAlgebra
using Plots
using StatsBase
using Printf
using ProgressBars
pythonplot()


# ## Noisy Single-Qubit Rabi Oscillations

# ### Solution to the master equation
# First, we observe the effect that incoherent depolarizing noise has on the Rabi oscillations of a single qubit. A depolarizing channel is modelled by collapse operators
# ``X, Y, `` and `` Z`` which occur with the a uniform rate ``\gamma``, as expressed in the following master equation
# ```math
# \rho(t+dt) = (1-idt\mathcal H)\rho(t)(1+idt H) + dt\gamma (X \rho(t) X + Y \rho(t) Y + Z \rho(t) Z-3\rho(t))
# ```
# Using the identity `` 2I = \rho + X\rho X + Y\rho Y + Z \rho Z`` for arbitrary `` \rho``, we can write
# ```math
# \rho(t+dt) = (1-idt\mathcal H') \rho(t) (1+idt \mathcal H') + 4\gamma dt\frac{I}{2}
# ```
# where ``\mathcal H' = \mathcal H - 4i\gamma``. This corresponds to ``\rho`` undergoing coherent evolution according to
# ``\mathcal H'`` with probability `` 1-4\gamma dt`` and being replaced with the completely mixed state 
# with probability `` 4\gamma dt``. Since all trace-preserving quantum channels stabilize the maximally mixed state, 
# the evolution can be modeled as a continuous-time Markov chain transitioning between coherent evolution
# ``|\psi(t)\rangle`` and the mixed state `` \frac{I}{2}`` with probability `` 4\gamma dt p_I(t)``, where `` p_I(t)`` 
# is the probability that the system is already in the mixed state. Integrating this probability over time gives ``p_I(t) = 1-e^{-4\gamma t}``. Therefore we can write down the solution
# ```math
# \rho(t) = e^{-4\gamma t}|\psi(t)\rangle\langle\psi(t)| + (1-e^{-4\gamma t})\frac{I}{2}
# ```math
# Where ``|\psi\rangle`` is evolved via the Schrodinger equation and normalized. Solving the Schrodinger equation for a 
# time-independent `` \mathcal H'`` gives `` |\psi(t)\rangle = e^{-i\mathcal Ht}e^{-4\gamma t}|\psi\rangle``. Normalizing
# this state gets rid of the exponential decay factor, leaving `` |\psi(t)\rangle\rangle = e^{-i\mathcal Ht}|\psi\rangle``, 
# corresponding to coherent evolution by the original Hamiltonian `` \mathcal H``.
#
# Consider the Hamiltonian ``\mathcal H = \frac{\Omega}{2}X`` and an initial state ``|\psi\rangle = |0\rangle``. 
# The time-evolution operator is ``e^{-iXt\Omega/2} = \cos(\frac{\Omega}{2}t) -i\sin(\frac{\Omega}{2}t)X``, 
# and so the state evovles in time to ``|\psi(t)\rangle = \cos(\frac{\Omega}{2}t)|0\rangle - i\sin(\frac{\Omega}{2}t)|1\rangle``.
# The expectation value of the number operator ``\hat n`` is then ``\langle \hat n(t) \rangle_\psi = \sin^2(\frac{\Omega}{2}t)``. 
# Plugging this into the master equation solution, in a noisy channel the evolution is
# ```math
# \langle \hat  n(t) \rangle = e^{-4\gamma t}\langle \hat n(t) \rangle_\psi + (1-e^{-4\gamma t})\frac{1}{2}\operatorname{Tr}(\hat n) = \frac{1}{2} - \frac{1}{2}e^{-4\gamma t}\cos(\Omega t)
# ```
# Like an underdamped harmonic oscillator, the value oscillates around an equilibrium value of ``1/2`` with an envelope that decays 
# exponentially in time.
#
# ### Simulation in Bloqade
# In BloqadeNoisy, the `NoisySchrodingerProblem` plays the same role as the `SchrodingerProblem` in Bloqade, and `emulate_noisy` plays the same role as `emulate!`. The `NoisySchrodingerProblem` is constructed using the following arguments:
# 1. initial state
# 2. list of times to save the solution
# 3. noiseless Hamiltonian
# 4. list of collapse operators. The rate is absorbed into the operators, so ``L_k`` becomes ``\sqrt{\gamma_k}L_k``. 
#
# `emulate` is called with the problem and the number of trajectories. Additional arguments can be a list of operators to take expectation values. Calling `emulate` with `report_error = true` will return a set of estimates of error by computing the standard deviation of the expectation values over the trajectories. Lastly, choices for the `ensemble_algo` argument are `EnsembleSerial()`, `EnsembleThreads()`, or `EnsembleDistributed()` for different levels of parallelization.

Ω = 2π
γ = 0.1

collapse_operator = sqrt(γ) .* mat.([X, Y, Z]) 
e_observable = [mat(Op.n)] #expectation values

h = rydberg_h([(0.0,)]; Ω = Ω)
save_times = LinRange(0, 10, 200)

prob = NoisySchrodingerProblem(
    zero_state(1), 
    save_times, 
    h, 
    collapse_operator
)

sol = emulate_noisy(
    prob, 2000, e_observable;
    ensemble_algo = EnsembleThreads(),
    report_error = true
)

plot(
    save_times,
    sol.expectations[1],
    ribbon = sol.twosigma[1],
    title = "Noisy Rabi Oscillation",
    xlabel = "Rabi periods",
    ylabel = L"\langle \hat n \rangle",
    label = "simulated"
)
plot!(
    save_times,
    1/2 .- 1/2*exp.(-4γ * save_times) .* cos.(Ω * save_times),
    label = "analytic"
)
ylims!(0,1)




# ![BloqadeNoisy](../../../assets/BloqadeNoisy_tutorial/noisy_depolarizingrabi.png)

# ## Coherent noise in neutral atom simulators
# One of the advantages of neutral atoms is that they couple weakly to their environment. This means that incoherent noise is supressed. A more dominant source of noise is due to imperfect control, which causes fluctuation the Rabi power ``\Omega``, detuning ``\Delta``, and atoms positions ``\vec r_{ij}`` between shots. This is referred to as "coherent" noise because it preserves the coherence over a single shot, but averaging over many shots still produces a mixed state. This noise afffects the system globally and is distinct from noise due to coupling to a bath.
#
# To model this affect, we consider ``\Omega`` distributed according to ``G(\Omega) = \frac{1}{\sqrt{2\pi\sigma^2}}e^{-(\Omega - \bar \Omega)^2/2\sigma^2}``. Let the Hamiltonian again be ``\mathcal H = \frac{\Omega}{2}X``. The noiseless evolution of ``\langle Z \rangle`` is then ``\langle Z(t) \rangle = \cos(\Omega t)``. We then average ``\langle Z \rangle`` over the distribution in ``\Omega``:
# ```math
# \overline{\langle Z \rangle} = \frac{1}{\sqrt{2\pi\sigma^2}}\int_{-\infty}^{\infty}\cos(\Omega t)e^{-(\Omega-\bar\Omega)^2/2\sigma^2} d\Omega = \cos(\bar \Omega t)e^{-2\sigma^2t^2/4}
# ```
# Thus we see where incoherent noise results in an exponential envelope, coherent noise results in a Gaussian envelope. This result holds generally. Combining this with the result from pure depolarizing noise, we obtain an envelope in the form of an exponentially modified Gaussian:
# ```math
# \langle Z(t)\rangle = \cos(\bar \Omega t)e^{-\sigma^2t^2/2 - 4\gamma t}
# ```
#
# ### Implementation in BloqadeNoisy
#
# The stochastic wavefunction method allows for this experimental setting to be modelled by varying the parameters in the Hamiltonian between trajectories. In BloqadeNoisy, this is represented as a method which takes a Hamiltonian and returns another method that can be called to generate random samples. The levels of coherent, incoherent, and readout noise are contained in an `ErrorModel` object, which is constructed using
# 1. A method that takes a number of atoms and returns a confusion matrix
# 2. A method that takes a number of atoms and returns a set of collapse operators
# 3. A method that accepts a Hamiltonian and returns a method to generate random samples

γ = .01
σ = .1 * 2π
Ω = 2π
h = rydberg_h([(0.0,)]; Ω = Ω)
times = LinRange(0,5,400)

em = ErrorModel(
    n->I,
    n->[sqrt(γ) * mat(op) for op in [X, Y, Z]], 
    h->(()->rydberg_h(
        [(0.0,)]; 
        Ω = Ω + σ*randn())
    ) 
)

prob = NoisySchrodingerProblem(zero_state(1), times, h, em)
sim = emulate_noisy(prob, 2000, [mat(Z)]; report_error = true);

plot(
    times, sim.expectations[1], ribbon = sim.twosigma[1], 
    label = "simulated", ylabel = L"\langle Z \rangle",
    xlabel = "Rabi periods", title = "Coherent noise"
)
plot!(
    times, 
    (t->cos(Ω*t)*exp.(-σ^2*t^2/2 - 4γ*t)).(times), 
    label = "analytic"
)

# ![BloqadeNoisy](../../../assets/BloqadeNoisy_tutorial/noisy_coherentnoise_rabi.png)

# ## Noise on Aquila
# Aquila's noise model was calibrated using experimental data, and the parameters that are modelled are the following:
#
# | Parameter  | Desc. | Default value |
# | ---------- | ----- | ------------- |
# | ``\delta \Omega/\langle \Omega \rangle`` | Relative global Rabi drive fluctation | 0.018 (dim.) |
# | ``\delta \Omega_i / \langle \Omega \rangle`` | Relative local Rabi drive fluctation | 0.008 (dim.) |
# | ``\delta \Delta`` | Global detuning fluctuation | ``0.18 \ (rad / \mu s)`` |
# | ``\delta \Delta_i`` | Global detuning fluctuation | ``0.50 \ (rad / \mu s)`` |
# | ``\delta r_x, \delta r_y`` | Lattice site position uncertainty | ``0.05 \ \mu m`` |
# | ``\tau_{relaxation}`` | Relaxation time | ``100 \ \mu s`` |
# | ``\tau_{dephasing}`` | Dephasing time | ``50 \ \mu s`` |
#
# In addition, the affect of readout noise can also be modelled using BloqadeNoisy. Readout error in neutral atom simulators is highly local, so it can be modeled as a problility ``p_{0 \to 1}`` and ``p_{1 \to 0}`` of confusing ``|0\rangle`` and ``|1 \rangle`` and vice-versa. 
#
#
# | Parameter  | Desc. | Default value |
# | ---------- | ----- | ------------- |
# | ``p_{0 \to 1}`` | Probability of confusing 0 with 1 | 0.01 |
# | ``p_{1 \to 0}`` | Probability of confusing 1 with 0 | 0.08 |
#
# The effect of readout error can be modelled by a confusion matrix of the form ``C = \begin{pmatrix} 1-p_{01} & p_{10} \\ p_{01} & 1-p_{10} \end{pmatrix}^{\otimes N}``. This transforms the output probability distribution from ``p(z)`` to ``\tilde p(z) = \sum_{x} C_{zx}p(x)``, which affects the sampled bitstrings and expectation values of operators estimated from the machine.
#
# ### Using Aquila noise model
# BloqadeNoisy makes Aquila's noise model available to use. It is constructed by calling `Aquila()` and is passed as the fourth argument to `emulate_noisy`.

# ### Experimental validation
# The experimental data used to calibrate the noise model is taken from the Aquila whitepaper [3]. Readout error can be added to expectation values in the computational basis by passing `readout = true` to `emulate`. The operators must be of type `Diagonal`. Below, we show the estimation of a noisy expectation value incorporating the effect of readout error and compare to the experimental whitepaper data.


dir = @__DIR__
filedir = joinpath(dir, "data", "15MHz_long.csv")
whitepaper_data = CSV.read(filedir, DataFrame, delim = ",", header = false)



times = collect(whitepaper_data[1,:])
data = collect(whitepaper_data[2,:])
save_times = LinRange(0, last(times), 400)

Ω = 15.3
Δ = 0
H = rydberg_h([(0.0,)]; Ω = Ω, Δ = Δ)

prob = NoisySchrodingerProblem(zero_state(1), save_times, H, Aquila())
sim = emulate_noisy(prob, 1000, [mat(Op.n)]; readout_error = true);

plot(times, data, marker = :diamond,
    markersize = 4, color = :green, linestyle = :dash,
    xlabel = L"t \ (\mu s)", ylabel = L"\langle \hat n \rangle",
    label = "Data", title = "Simulation vs Data"
)
plot!(save_times, sim[1], color = :blue, label = "Sim")

# ![BloqadeNoisy](../../../assets/BloqadeNoisy_tutorial/noisy_whitepapercomparison.png)

# The `Aquila()` method is shorthand for `load_error_model(JSON.parse(AQUILA))`. The `AQUILA` string is a JSON dictionary containing the preconfigured noise model, and this syntax allows the noise model to be modified.

default_error_model = JSON.parse(AQUILA)
dephasing_rate = 1/10
default_error_model["incoherent"]["dephasing rate"] = dephasing_rate

#create error model based on Aquila
new_error_model = load_error_model(default_error_model);

# ### Simulated scar experiment
#
# Many-body scars are product states that have a large overlap with a tower of eigenstates exhibiting atypical properties such as low entanglement. Ref [4] studies the scarring properties of the PXP model, given by the Hamiltonian
# ```math
# \mathcal H_{\text{PXP}} = \sum_{i}P_iXP_i + \mu \sum_{i}\hat n_i
# ```
# where ``P_i = \bigotimes_{\langle j, i \rangle}|0_i\rangle\langle 0_i|``, with ``\langle j, i\rangle`` running over the nearest neighbors of ``i``, and ``\mu`` is the chemical potential. Scarring is observed when states fail to thermalize, i.e. exectation values exhibit long-lived dynamic revivals. The PXP model can be approximately realized in a Rydberg atom array by operating in the blockaded subspace and tuning ``\Delta`` to cancel the mean-field contribution of the next-nearest-neighbor Van der Waals interactions. 
#
# BloqadeNoisy allows for estimation of expectation values at a fixed number of shots to emulate an experimental setting. We can use this to predict what a many-body scar experiment on Aquila would look like. Following [4], we first adiabatically prepare the ground state of the PXP model with a chemical potential ``\mu_i = -0.76`` and then a quench to ``\mu = 1.6``. 

Ω = 2π
C = 862690 * 2π
Rb = (C/Ω)^(1/6)
N = 8

atoms = generate_sites(ChainLattice(), N, scale = .65*Rb)
Δ_NNN = 1/(2N) * sum(
    [ 
        (r=norm(atoms[i] .- atoms[j]);
        r > Rb ? C/r^6 : 0)
        for i in 1:N for j in 1:N
    ]
)

μ_i = -.76 * 2π
μ_f = 1.6 * 2π

t_prep = 2 
t_quench = 2 

Δ_ad = piecewise_linear(clocks = [0, t_prep], values = [-40*2π, Δ_NNN - μ_i])
Δ_quench = piecewise_constant(clocks = [0, t_quench], values = [Δ_NNN - μ_f])
Δ = append(Δ_ad, Δ_quench)

times = LinRange(t_prep, t_prep+t_quench, 200)
H = rydberg_h(atoms; Ω = Ω, Δ = Δ)
ψ = solve(
    SchrodingerProblem(zero_state(N), t_prep+t_quench, H;
        save_on = true, saveat = times, save_start = false
    ),
    DP8()
);

expt_times = LinRange(t_prep, t_prep+t_quench, 70)
prob = NoisySchrodingerProblem(zero_state(N), expt_times, H, Aquila())
@time sim = emulate_noisy(prob, 200, [mat(chain(igate(N)-put(N, i=>Op.n) for i in 1:N))];
    readout_error = true, report_error = true, shots = 1000,
    ensemble_algo = EnsembleThreads()
)

plot(times, [abs(p[1])^2 for p in ψ.u],
    label = "ideal",
    color = :blue, title = "Many-body scar simulation"
)
scatter!(expt_times, sim.expectations, 
    yerr = sqrt.(sim.shot_error[1].^2 .+ sim.propagated_err[1].^2),
    xlabel = L"t \ (\mu s)",
    ylabel = L"\langle (1-\hat n)^{\otimes N}\rangle",
    label = "noisy", color = :red,
)

# ![BloqadeNoisy](../../../assets/BloqadeNoisy_tutorial/noisy_manybodyscar.png)




# ## References & Further Reading
# [1] https://qutip.org/docs/latest/guide/dynamics/dynamics-monte.html
#
# [2] https://lukin.physics.harvard.edu/files/lukin/files/physics_285b_lecture_notes.pdf (Chapter 6)
#
# [3] Wurtz, J. et al. (2023). Aquila: QuEra's 256-qubit neutral-atom quantum computer. arXiv preprint arXiv:2306.11727.
#
# [4] Daniel, Aiden, et al. "Bridging quantum criticality via many-body scarring." Physical Review B 107.23 (2023): 235108.
