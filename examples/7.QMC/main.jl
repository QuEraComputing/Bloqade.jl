# # Quantum Monte Carlo Method
# 
# You remember that tutorial about adiabatic evolution? The one in which we prepared ordered ground states of Rydberg Hamiltonians, such as the $\mathbb{Z}_2$  phase? What method did we use to achieve that goal? 
# Not sure? Don't worry - we actually never mentioned its name. That's because that tutorial was all about introducing those interesting ground states themselves. In this tutorial, on the other hand, the aim is to familiarize ourselves with a specific *method* that can produce such results - the **Quantum Monte Carlo** method. In particular, it will allow us to tackle system sizes far beyond what is possible with exact diagonalization - yup, that's the method used in the previous tutorial. While ED is limited to a few tens of particles, QMC can deal with hundreds.
# Before we dive into QMC, let's take a step back and make sure we understand the gist of plain Monte Carlo. Consider the following example: You are asked to calculate a really hard integral. You have no clue how to solve it analytically. Is there another way you could approach the problem? 
# You could draw the corresponding graph, add a square around it and throw darts randomly at the square. For each throw, you record whether the dart landed below or above the graph. Then you form the ratio $N_{below}/N_{above}$ and multiply it by the area of the square. Indeed, after enough throws that value will approach the true area underneath the graph.
# 
# ![Integral](../../../assets/QMC_tutorial/integral.png)
# 
# *Note to myself: Add in a square to image*

# Of course, you might spare yourself the darts and ask a computer to generate a very large number of uniformly random points scattered across the plot, to improve the accuracy of your result. As a matter of fact, the history of the Monte Carlo method and the rise of computers are closely intertwined - the very [first MC](https://en.wikipedia.org/wiki/Monte_Carlo_method#History) was run on the ENIAC itself in 1948 to simulate neutron diffusion processes in the hydrogen bomb. It was this MC simulation that gave birth to the modern era of computational physics. 
# To summarize, what's the gist of Monte Carlo? Instead of solving a problem exactly, you invoke randomness to sample from the distribution involved in the problem, in this case the $f(x)$ in $\int_0^3 f(x) dx$, until your result approximates the true solution closely enough. This begs the question of what defines *closely enough*, an issue we will examine in more detail below.

# ### Quantum Monte Carlo
# 
# So what about *Quantum* Monte Carlo? Firstly, let us emphasize that QMC is still a *classical* simulation, i.e. it is not run on a quantum machine such as QuEra's Aquila. It is simply called *quantum* because we aim to use the idea of MC to investigate problems from quantum physics. In other words, the space we are generating samples from is a Hilbert space of quantum-mechanical configurations. Furthermore, QMC is one of the best established methods in numerically tackling the analytically intractable integrals of quantum *many-body* physics that are beyond the reach of exact solutions. 
# Typically, the integral or sum in question is the expection value of some observable $\langle A \rangle_\psi = \sum_j a_j |\langle \psi | \phi_j \rangle |^2$ such as the energy, magnetization etc. The issue we face is the *curse of dimensionality*, meaning the number of terms in this sum grows exponentially in system size. We circumvent this curse by not calculating the whole sum, but instead sampling from the probability distribution given by $|\langle \psi | \phi_j \rangle |^2$, favoring those terms that contribute signficantly to the sum, i.e. have large weights $a_j$. While doing so, it is essential that we explore the configuration space ergodically. This simply means that configurations with a small but non-zero weight should have a chance of being reached, even though this will occur less frequently than for those with large weights. 
# To finish of our introduction, we should mention that QMC is an umbrella term comprising many different implementations of this same core idea, each tailored to a specific class of quantum problems. Under the hood, BloqadeQMC currently implements the SSE (Stochastic Series Expansion) method. It was [first invented](http://physics.bu.edu/~sandvik/research/ssehistory.html) by Anders Sandvik to study e.g. Heisenberg-type models and recently adapted to the Rydberg Hamiltonian by [Merali et al (2021)](https://arxiv.org/abs/2107.00766).

# ### Getting Started with BloqadeQMC
# 
# Let's get started in the usual way by importing the required libraries:

using BloqadeQMC
using Random
using Plots

# In addition, we'll import some libraries that we can later use to check our QMC against ED results for small system sizes. 

using Bloqade 
using Yao: mat, ArrayReg
using LinearAlgebra
using Measurements
using Measurements: value, uncertainty
using Statistics


# As an initial example, let's recreate the $\mathbb{Z}_2$ phase in the 1D chain which we first saw in the tutorial on *Adiabatic Evolution*. This means defining the lattice geometry as well as the Rabi drive and global detuning parameters which we'll once again feed into the Rydberg Hamiltonian. 

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.72)

Ω = 2π * 4
Δ = 2π * 10
h = rydberg_h(atoms; Δ, Ω)


# Now comes the first new step:

h_qmc = rydberg_qmc(h)

# What has happened in this step? The object h_qmc still contains all the previous information about the lattice geometry as well as the Hamiltonian parameters $\Omega$ and $\Delta$. Crucially, however, this object now also stores the distribution of weights from which our algorithm will sample. Without going into all the details which can be found in [Merali et al (2021)](https://arxiv.org/abs/2107.00766), let's focus on those key elements of the SSE formalism that you will need to calculate observables from the samples. This means understanding what exactly is meant by configuration space in the SSE formalism, what samples from that space look like and what their weights are. 

# ### Prelude: Massaging the partition function
# Before answering those questions, let us revisit the finite temperature partition function $Z = Tr(e^{-\beta H})$. Indeed, $Z$ is the protagonist in the mathematical formalism of SSE. Massaging it through a few tricks and combinatorics will help us answer the questions and prepare us for the picture that will come below.
#
# The core idea is the following: Instead of calculating the trace analytically, we can first write out the Taylor series of the exponential. (That's where the SE in SSE, the idea of *series expansion*, comes in!).
# $$
#     Z = Tr(e^{-\beta H}) = Tr(\sum_{n=0}^{\infty} \frac{\beta^n}{n!}(-\hat{H})^n).
# $$
# Next, we simply insert the usual identities over a set of basis states $\{\alpha_p\}$, giving
# $$
# Z = \sum_{\{\alpha_p\}} \sum_{n=0}^{\infty}\frac{\beta^n}{n!} \prod_{p=1}^n \langle\alpha_{p-1}| -\hat{H} | \alpha_p\rangle.
# $$
# Finally, we split the Hamiltonian into a sum of local Hamiltonians where each term acts on only one or two atoms. (For example, the latter includes the Rydberg interaction term.) Formally, we'll write $H = - \sum_{t, a} \hat{H}_{t,a}$  where the labels $t$ symbolizes whether the term is diagonal ($t=1$) or off-diagonal ($t=-1$) and $a$ specifies the atoms the term acts on. This leads to the last crucial step in the massage: We switch the order of the n-fold product and sum. From each copy of the total Hamiltonian, we pick one of the local terms and multiply these in a new n-fold product:
# $$
#     Z =  \sum_{\{\alpha_p\}} \sum_{n=0}^{\infty} \sum_{S_n} \frac{\beta^n}{n!} \prod_{p=1}^n \langle\alpha_{p-1}| - \hat{H}_{t_p, a_p} | \alpha_p\rangle.
# $$
# *Note: For further technical details, please refer to section 2.1 of [Merali et al (2021)](https://arxiv.org/abs/2107.00766)*
# Why is this crucial? This procedure of picking local terms produced something we will henceforth refer to as an *operator sequence* and denote by $S_n$ with $n$ being the length of this sequence. Let's have a look what this looks like!

# ## The SSE Configuration Space
# Let's consider a system consisting of four Rydberg atoms. One state in the SSE configuration space of this system could look like the following,

#![SSE](../../../assets/QMC_tutorial/SSEconfig.png)
# where we will briefly define each ingredient in this probably surprising picture. 
#
# One aspect of the picture we can directly identify is its vertical dimension: the four circles correspond to four atoms, with filled circles being excited atoms $|r\rangle$ and white circles representing the ground state $|g\rangle$.
# 
# Next, let's define what the blue, red and black boxes situated either on one horizontal line or sitting between two lines stand for. These are the local Hamiltonian terms we mentioned earlier. The interested reader may find the full definition of the local terms in section 3 of [Merali et al (2021)](https://arxiv.org/abs/2107.00766).
# 
# Next, let's try to understand the horizontal dimension of the picture: Why are there seven copies of the four atoms? Or equivalently, why are there 6 boxes? This is the visual equivalent of looking at the sixth expansion order in the Taylor series introduced above! Indeed, each SSE *sample* will correspond to a particular expansion order $m$ since the infinite Taylor series gave rise to an intractable sum which we can use the core idea of MC to sample from. Neat, right? It is that number $m$ that defines the horizontal length of the picture.
# 
# Putting the pieces together, we are ready to get the gist of the weight distribution stored in the rydberg_qmc object. That distribution is given by the matrix elements defined by the local operators (boxes) sandwiched in between two configurations of the physical system (two copies of the four vertically stacked atoms). The picture as a whole is what defines one sample in the SSE configuration whose overall weight is given by the product of the matrix elements of each operator.
# 
# Don't worry if it takes you a while to process the above information. In fact, a precise understanding is not necessary to successfully run a simulation using BloqadeQMC. We have simply included this peek into the backend in order for you to have some notion of what we mean when we refer to the *number of operators sampled* in the energy calculations later on. We simply mean the number of boxes in Figure xxx.

# ### Running a simulation using BloqadeQMC
# 
# Now that we have defined the configuration space, let's traverse it and generate samples from it. Importantly, there is a key difference in how these samples are generated as opposed to when we were throwing darts in the beginning. Those throws were all individually random and independent. The algorithm implemented in BloqadeQMC, on the other hand, falls within the class of Markov Chain Monte Carlo (MCMC) methods. As the name suggests, the samples are no longer independent but form a chain in which the probability of the next sample depends on the current sample. (This is the so-called Markov property.)

# Before jumping into the code, let's touch on the prototypical example of how such a chain might be built: via the Metropolis-Hastings algorithm. This recipe has three main steps:
# 1) randomly propose a new configuration, 
# 
# 2) accept the proposal according to a probability given by the ratio of the weights of the current and proposed configuration, 
# $$P = min(\frac{W_{current}}{W_{proposed}}, 1)$$
# 
# 3) repeat.

# Now, let's define the parameters than govern the length of the chain, i.e. how long we let the simulation run for.

EQ_MCS = 100
MCS = 100_000

# As you can see, we are in fact not defining one but two parameters. Indeed, a MCMC simulation is typically split into two phases: first the equilibration phase, also referred to as the *burn-in*, followed by the sampling phase. Even though the mathematical theorems surrounding Markov chains guarantee that under reasonable assumptions these chains eventually converge to the desired probability distribution, it does take some time until this is indeed the case.
# 
# While there is no formula that will tell you when precisely equilibration has been reached, a good heuristic is plotting an observable such as the energy. You will see its value initially fluctuating and eventually stabilizing around a value which you can then assume is its equilibrium expectation value.*
# 
# ![Equilibration](../../../assets/QMC_tutorial/equilibration.png)
# Now, we're almost ready to run the simulation. You can think of the BinaryThermalState object as the initialization (rephrase). You can roughly think of M as the maximum expansion order. (Indeed, we are allowed to truncate the expansion since one can show that the weights of higher terms in the expansion fall off as a Poisson distribution.) 
# 
# We choose the inverse temperature β large enough for the simulation to approximate the ground state. If your system is gapped, then a finite value for β will be enough to reach the ground state. If the gap closes, then you will need to scale β with the system size.
M = 50
ts = BinaryThermalState(h_qmc, M)   
# BinaryThermalState is an object necessary in the backend to store the instantaneous SSE configuration during the MC steps.
d = Diagnostics()                   
# Diagnostics are a feature that can be used by the advanced user to analyse performance and extract further information from the backend. We refer to the manual (in progress) for details.


rng = MersenneTwister(3214)

β = 0.5


# Time to run the simulation!
#
# *Note: You will see that mc_step_beta!() returns three objects which can be used for computations during each MC step. "ts" stores the instantaneous SSE configuration, "h_qmc" is the same object as before, "lsize" is related to the arrays used to the arrays used to carry out a MC step and depend on the operator sequence of the current configuration. Details are explained in the upcoming manual. Further, SSE_slice stores a sample of the atom configuration taken from the current SSE configuration, in this case chosen to be the first vertical slice."

[mc_step_beta!(rng, ts, h_qmc,β, d, eq=true) for i in 1:EQ_MCS] # equilibration phase

densities_QMC = zeros(nsites)
occs = zeros(MCS, nsites)

for i in 1:MCS # Monte Carlo Steps
    mc_step_beta!(rng, ts, h_qmc,β, d, eq=false) do lsize, ts, h_qmc
        SSE_slice = sample(h_qmc,ts, 1)
        occs[i, :] = ifelse.(SSE_slice .== true, 1.0, 0.0)
    end
end

for jj in 1:nsites
    densities_QMC[jj] = mean(occs[:,jj])
end


# Let's plot the results.

using Plots: bar
    
bar(densities_QMC, label="")
xlabel!("Site number")
ylabel!("Occupation density")


# <img alt="9atoms" src="./9atom_chain_Z2.png">
# ![9atoms](../../../assets/QMC_tutorial/9atom_chain_Z2.png)

# As expected, we see a $\mathbb{Z}_2$ pattern has emerged, just as we saw using the exact diagonalization method. So let's try an example that goes beyond what is feasible with ED. We run the same code as before, substituting the 1D chain with 9 atoms for a 2D square lattice with 100 atoms. (This should only take a minute or two to run on your laptop.)

nx = ny = 10
nsites = nx*ny
atoms = generate_sites(SquareLattice(), nx, ny, scale = 6.51);

# ![100atoms](../../../assets/QMC_tutorial/Checkerboard_10x10.png)

# ### Calculating observables using BloqadeQMC
# 
# The final question we will address in this tutorial is how to calculate observables using BloqadeQMC. We'll choose to investigate the energy during a detuning sweep. Furthermore, we'll limit ourselves to a system size which ED can also handle such that we may compare the results from both methods.
# 
# Let's start by returning to the 1D chain and defining the Hamiltonian parameters. We will keep the constant Rabi drive from before, i.e. $\Omega = 4 \times 2\pi$ MHz. For the detuning, we will choose the following ramp.

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.72);

Δ_step = 15
Δ = LinRange(-2π * 9, 2π * 9, Δ_step)


# Now, we only need to make two small changes to our previous MC code. First, for each value of $\Delta$ specified in this ramp, we will run a separate QMC simulation that will produce one data point in the final energy plot. Second, we need to store the number of operators contained in each sample. Yes, this is where the picture we introduced earlier comes in. This number will fluctuate from sample to sample as we build the Markov chain. It is directly returned by the mc_step_beta!() function.
# 
# We can then use the number of operators to compute the energy thanks to the following observation. Firstly, we know from statistical physics that
# 
# $$
#     \langle E \rangle = -\frac{\partial ln Z}{\partial \beta}.
# $$
# 
# Yet in the SSE formalism we can find another expression for the expectation value $\langle E \rangle$. Let's work backwards and consider the SSE expectation value of the length of the operator sequence:
# 
# $$
#     \langle n \rangle = \frac{1}{Z} Tr(\sum_{n=0}^\infty n \frac{(-\beta \hat{H})^n}{n!}) \\
#     = \frac{1}{Z} Tr(\sum_{n=1}^\infty n \frac{(-\beta \hat{H})^n}{n!}) \\
#     = \frac{1}{Z} Tr(\sum_{n=1}^\infty \frac{(-\beta \hat{H})^n}{(n-1)!}) \\
#     = \frac{1}{Z} Tr(\sum_{n=0}^\infty (-\beta \hat{H}) * \frac{(-\beta \hat{H})^n}{n!}),
# $$
# 
# where in the second line we dropped the vanishing $n=0$ term and in the fourth line we shifted the sum over n by 1 and hence obtained the extra factor $-\beta \hat{H}$. From this, we directly read of the SSE formula for the energy expectation value:
# 
# $$
#     \langle E \rangle = -\frac{\langle n \rangle}{\beta}.
# $$

using BinningAnalysis

energy_QMC = []

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
    h_ii_qmc = rydberg_qmc(h_ii)
    ts_ii = BinaryThermalState(h_ii_qmc,M)
    d_ii = Diagnostics()

    [mc_step_beta!(rng, ts_ii, h_ii_qmc, β, d_ii, eq=true) for i in 1:EQ_MCS] #equilibration phase
    
    ns = zeros(MCS)

    for i in 1:MCS # Monte Carlo Steps
        ns[i] = mc_step_beta!(rng, ts_ii, h_ii_qmc,β, d_ii, eq=false)
    end

    energy(x) = -x / β + h_ii_qmc.energy_shift  # The energy shift here ensures that all matrix elements are non-negative. See Merali et al. for details.
    BE = LogBinner(energy.(ns)) # Binning analysis 
    τ_energy = tau(BE)
    ratio = 2 * τ_energy + 1
    energy_binned = measurement(mean(BE), std_error(BE)*sqrt(ratio)) 
    append!(energy_QMC, energy_binned)
end


# *Note: For a detailed guide to the backend of the mc_step_beta!() function, please see the manual (in progress).*

# We must address one last point when calculating expectation values using a MCMC method. As discussed, the samples are not statistically independent, instead they are correlated. This point is crucial. After all, Monte Carlo is an approximate method. We must construct error bars around the mean values generated by our simulation. Correlation of samples has the effect of reducing the effective variance, thus purporting an accuracy we cannot in fact justify. 
# 
# One method of addressing this issue is via binning analysis. The idea is to bin small sequences of samples together and average over them. This gives rise to an reduced effective sample size given by $N_{eff} = \frac{N_{orig}}{2\tau + 1}$, where $\tau$ is the autocorrelation time. We use this effective sample size to rescale the standard error by $\sqrt{2\tau + 1}$. For details, please see the documentation of the imported BinningAnalysis package. If the error bars thus achieved are still too large, try increasing the number of samples taken, i.e. the number of Monte Carlo steps.

# Finally, let's carry out the ED calculation as well and plot both results together.

energy_ED = zeros(Δ_step)

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
    h_m = Matrix(mat(h_ii)) 
    energies, vecs = LinearAlgebra.eigen(h_m) 

    w = exp.(-β .* (energies .- energies[1]))
    energy_ED[ii] = sum(w .* energies) / sum(w)
end

scatter(Δ/2π, energy_ED, label="ED", marker=:x)
scatter!(Δ/2π, value.(energy_QMC); yerror=uncertainty.(energy_QMC), label="QMC", marker=:x)
xlabel!("Δ/2π")
ylabel!("Energy")

# ![energy](../../../assets/QMC_tutorial/Energy_Comparison.png)

# We see that using the QMC, we have indeed achieved the same results as for the ED with high accuracy.
# 
# *Note: It is possible to also calculate other observables in the SSE formalism, if they are products of the operators present in the Hamiltonian itself. For example, to calculate $\langle X \rangle$, one would count the number of times the off-diagonal operator $\sigma_x$ appears in each sample and take its average. For more examples and derivations, please see [insert Sandvik 90s paper]*

# To finish off this tutorial, we will leave with one final plot of the staggered magnetization which acts as an order parameter to observe the transition from the disordered to the $\mathbb{Z}_2$ phase achieved before :)

densities_QMC = []
order_param_QMC = []

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
    h_ii_qmc = rydberg_qmc(h_ii)
    ts_ii = BinaryThermalState(h_ii_qmc,M)
    d_ii = Diagnostics()

    [mc_step_beta!(rng, ts_ii, h_ii_qmc, β, d_ii, eq=true) for i in 1:EQ_MCS] #equilibration phase
    
    order_param = zeros(MCS)

    for i in 1:MCS # Monte Carlo Steps
        mc_step_beta!(rng, ts_ii, h_ii_qmc, β, d_ii, eq=false) do lsize, ts, h_qmc
            SSE_slice = sample(h_qmc,ts, 1) # occ = 0,1
            spin = 2 .* SSE_slice .- 1 # spin = -1, 1
            order_param[i] = abs(sum(spin[1:2:end]) - sum(spin[2:2:end]))/length(spin)
        end
    end

    BD = LogBinner(order_param)
    τ_energy = tau(BD)
    ratio = 2 * τ_energy + 1
    energy_binned = measurement(mean(BD), std_error(BD)*sqrt(ratio)) 
    append!(order_param_QMC, measurement(mean(BD), std_error(BD)*sqrt(ratio)) )
end

scatter(Δ/2π, value.(order_param_QMC); yerror=uncertainty.(order_param_QMC), label="", marker=:x)
xlabel!("Δ/2π (MHz)")
ylabel!("Stag mag")

# ![order_param](../../../assets/QMC_tutorial/order_parameter.png)
