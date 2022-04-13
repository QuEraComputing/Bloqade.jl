# # Background
# In [H. Pichler, et al.](https://arxiv.org/pdf/1808.10816.pdf), it was shown that 
# Rydberg atom arrays can be used to encode the maximum independent set (MIS)
# problem on unit disk graphs (UDG).  In this example, we present an implementation of 
# quantum annealing to solve the maximum weight independent set (MWIS) problem 
# of a weighted unit disk graph, with arbitrary weights for each vertex.  The MWIS problem 
# seeks to find the independent set whose weights sum to the maximum possible value. 

# We import the required packages to compute weighted MIS classically

using Graphs, GraphTensorNetworks

# We initially specify the atom locations and construct the corresponding diagonally-coupled 
# unit disk graph on a square lattice.  The atoms represent vertices on the problem graph, 
# and all vertices closer than a length 1.5 are connected by an edge.  
locs = [(1,-1), (4,0), (1,1), (2,0), (0,0), (2,2), (2,-2), (3,1), (3,-1)];
g = unit_disk_graph(locs, 1.5)
show_graph(g; locs=locs, vertex_colors=["white" for i=1:nv(g)])

# We assign random weights to each vertex.  
weights = [rand() for i = 1:nv(g)];

# We solve the MWIS problem classically for this graph for reference 
# using [GraphTensorNetworks package](https://github.com/QuEraComputing/GenericTensorNetworks.jl). 
# The MWIS is shown in red.
configs_mapped = solve(IndependentSet(g; weights= weights), ConfigsMax())[]
MIS_config = configs_mapped.c[1]
show_graph(g; locs = locs, vertex_colors=
          [iszero(MIS_config[i]) ? "white" : "red" for i=1:nv(g)])


# # Quantum Annealing to solve the MWIS problem

# A quantum annealing algorithm (QAA) can be performed with the Hamiltonian:

# $H_{QA}(t) = \sum_{v \in V} (- \Delta_v(t) n_v + \Omega_v(t) \sigma_v^x) + \sum_{(u, w) \in E} U_{u,w} n_u n_w$

# We will work in the limit of $\Delta, \Omega \ll U$, where 
# the non-independent set space of the graph can be neglected (in experiments, this corresponds 
# to the limit where the interaction energy is much stronger than other energy scales).  
# In this limit, we can restrict outselves to the Rydberg blockade subspace (see [blockade](@ref)) of the graph 
# and ignore the interaction term $\sum_{(u, w) \in E} U n_u n_w$ in the Hamiltonian. We can 
# vertex weights of the MWIS problem in this Hamiltonian with individual atom detuning (specifying $\Delta_v(t)$
# for each atom)


# The QAA can be designed by first initializing all qubits to the ground state of 
# $H_{QA}(t = 0)$ where $\Delta(t = 0) = 0$ and $\Omega(t = 0) = \Omega_{max}$.
#  We then change the parameters by turning down $\Omega(t)$ to 0, and turning 
#  up $\Delta(t)$ to $\Delta_{max}$ after some final time $t_f$.  

# By the adiabatic theorem, when the time evolution is sufficiently slow, 
# the system follows the instantaneous ground state and ends up in the 
# solution to the MWIS problem.  


# ## Build Pulse Sequence
# Because we are considering the weighted MIS problem, we implement individual
# atom detuning with $\Delta(t)_i = w_i \times \Delta(t)$.  We can simulate individual 
# pulse shapes on the emulator.

using Bloqade, BloqadePlots
using PythonCall
plt = pyimport("matplotlib.pyplot")

# We build the Hamiltonian and the corresponding waveforms for adiabatic evolution of the system
function build_adiabatic_sweep(graph, Ω_max::Float64, Δ_max::Float64, t_max::Float64, weights::Vector{Float64})
    Ω = Waveform(t->Ω_max + (- Ω_max)*t/t_max, duration=t_max)
    Δ = map(1:nv(graph)) do idx
        Waveform(t->weights[idx] * Δ_max * t/t_max, duration=t_max)
    end

    h = XTerm(nv(graph), Ω) - NTerm(nv(graph), Δ)
    return h, Ω, Δ
end;

# We choose $\Delta_{max} / \Omega_{max} = 2.5$, with $\Omega_{max} = 4 \times 2 \pi$ MHz
Ω_max = 2 * 4 * pi
Δ_max = 2.5 * Ω_max
t_max = 1.0
h, Ω, Δ = build_adiabatic_sweep(g, Ω_max, Δ_max, t, weights);

fig, (ax1, ax2) = plt.subplots(nrows=2)
draw!(ax1, Ω) 
for i = 1:nv(g)
    draw!(ax2, Δ[i])
end 
fig


# # Compute MIS Probability and Adiabatic Timescale

# We import additional libraries to solve for the ground state of the initial Hamiltonian
using SparseArrays

# We compute the MIS probability of the original graph as a function of time.  We 
# want to extract the adiabatic timescale $T_{LZ}$ from the Landau-Zener fitting: 
# $1 - P_{MIS} = e^{a - T/T_{LZ}}$.  To do this, we find the first instance time $T^*$
# such that P_{MIS}(T) > 0.9, and continue to run evolutions to 2.5T^* to extract $T_L{LZ}$.


# We work in the blockade (independent set) subspace

t_list = []
P_MIS = []
subspace = independent_set_subspace(g)

global t = 0.1
global T_  = 0.0

while (t < T_  * 2.5) || (T_  == 0.0)
    global t 
    global T_

    h = build_adiabatic_sweep(g, Ω_max, Δ_max, t, weights)[1]
    
    # compute the ground state of the initial Hamiltonian
    energies, GS = eigsolve(SparseMatrixCSC(h(0.0), subspace), 1, :SR);
    r = RydbergReg(GS[1], subspace);
    
    # run ODE evolution
    prob = ODEEvolution(r, t, h; dt=1e-3, adaptive=false)
    emulate!(prob)
    
    # compute MIS probability
    p = maximum_independent_set_probability(prob.reg, g, MIS_config)
    
    # find first occurrence p > 0.9
    ((p > 0.9) && (T_ == 0)) && (global T_ = t)
    
    push!(t_list, t)
    push!(P_MIS, p)
    
    if (p > 0.9)
        global t += T_ * 0.25
    else
        global t += 0.1
    end
end

# We can compute the adiabatic timescale by fitting a Landau Zener curve to the 
# MIS probability: 
using CurveFit
y = broadcast(log, 1 .- P_MIS[P_MIS .> 0.9])
a, b = linear_fit(t_list[P_MIS .> 0.9], y)
T_LZ = -1/b

# Plot results
fig, (ax1, ax2) = plt.subplots(ncols = 2)
ax1.scatter(t_list, 1 .- P_MIS)
ax1.set_ylabel("MWIS Probability")
ax1.set_xlabel("Time (μs)")

ax2.scatter(t_list, broadcast(log, 1 .- P_MIS))
ax2.plot(t_list, a .+ b .* t_list)
ax2.set_yscale('log')
ax2.set_xlabel("Time (μs)")
ax2.set_ylabel("1 - P(MWIS)")

fig
