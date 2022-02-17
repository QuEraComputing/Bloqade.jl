# # Background
# In [H. Pichler, et al.](https://arxiv.org/pdf/1808.10816.pdf), it was shown that 
# Rydberg atom arrays can be used to encode the maximum independent set (MIS)
# problem on unit disk graphs (UDG).  In this example, we present an implementation of 
# quantum annealing to solve the weighted maximum independent set (MIS) problem 
# beyond UDGs, for the specific example of the non-UDG graph: $K_{2,4}$.
# We use the emulator to analyze the performance of quantum optimization algorithms
# for the mapped UDG.  

# # Transform Arbitrary Graphs to Unit Disk Graphs 
# Using a series of encoding tools and at most a quadratic overhead, 
# we can map the weighted MIS problem of an arbitrary graph to the 
# weighted MIS problem of a corresponding generated, transformed UDG 
# on a square lattice. The MIS solution of the mapped graph will contain 
# the MIS solution of the original graph.

# We import the required packages for mapping the original graph and visualization tools
using Graphs, GraphTensorNetworks, UnitDiskMapping
using Plots, Plots.PlotMeasures, LaTeXStrings, Statistics
using SparseArrays, ThreadsX, BitBasis
import GraphTensorNetworks.visualize

# The original graph is a non-UDG. 
g_0 = SimpleGraph(6)
for (i,j) in [(5, 1), (5,2), (5,3), (5,4), (6, 1), (6, 2), (6, 3), (6, 4)]
    add_edge!(g_0, i, j)
end
pos_0 =  [(0,0), (1,0), (2,0), (3,0), (1, 1), (2,1)]
show_graph(g_0; locs = pos_0)

# We can use a series of gadgets to map the graph into a UDG 
M = map_graph(Weighted(), g_0, vertex_order=[6,5,4,3,2,1]);
g_udg = M.grid_graph

# We can apply simplification rules (will automate into UnitDiskMapping soon) to 
# simplify the mapped graph with coordinates at new_pos.  This generates a unit disk
# graph on the square lattice, where nodes closer than length 1.5 are connected 
# by an edge.  The nodes corresponding to the nodes in the original graph, where the MIS 
# of the mapped graph can be directly mapped back to that of the original graph, are 
# shown in blue. The white nodes represent ancilla vertices. 
simplified_locs = [ (6,1), (4,0), (2,0), (0,0), (1,1), (1,-1), 
 (2,2), (2,-2), (6, -1), (7,0), (3,1), (3,-2), (4,2), (4,-2), (5,1), (5,-1)]



g_m = GraphTensorNetworks.unit_disk_graph(simplified_locs, 1.5)
show_graph(g_m; locs=simplified_locs, vertex_colors=
          [i > nv(g_0)  ? "white" : "lightblue" for i=1:nv(g_m)])

# We solve the weighted MIS problem for a randomly weighted $K_{2,4}$ graph, normalized 
# such that the maximum weight is 1.0.  The corresponding mapped graph has the same weights 
# as the original graph, with the ancilla vertices all having weight 1.0. We 
# generate random weights for the original graph (between 0.5 and 1.0):
weight_Δ_M = ones(nv(g_m))
weight_Δ_0 = [rand() * 0.5 + 0.5 for i = 1:nv(g_0)]
max_weight = maximum(weight_Δ_0)

for i = 1:nv(g_0)
    weight_Δ_0[i] = weight_Δ_0[i] / max_weight
    weight_Δ_M[i] = weight_Δ_0[i]
end 

# # Quantum Annealing for $K_{2,4}$ graph

# A quantum annealing algorithm (QAA) for MIS can be performed with the Hamiltonian:

# $H_{QA}(t) = \sum_{v \in V} (- \Delta(t) n_v + \Omega(t) \sigma_v^x) + \sum_{(u, w) \in E} U n_u n_w$

# Because the original graph is not a UDG, and thus do not have a natural encoding in the 
# Rydberg hardware geometry, we will work in the limit of $\Delta, \Omega \ll U$, where 
# the non-independent sets of the graph can be neglected (in experiments, this corresponds 
# to the limit where the interaction energy is much stronger than other energy scales).  
# Thus, in this limit, we can restrict outselves to the Rydberg blockade subspace of the graph 
# and ignore the interaction term $\sum_{(u, w) \in E} U n_u n_w$ in the Hamiltonian. 

# The QAA can be designed by first initializing all qubits to the ground state of 
# $H_{QA}(t = 0)$ where $\Delta(t = 0) = 0$ and $\Omega(t = 0) = \Omega_{max}$. 
#  We then change the parameters by turning down $\Omega(t)$ to 0, and turning 
#  up $\Delta(t)$ to $\Delta_{max}$ after some final time $t_f$.  

# By the adiabatic theorem, when the time evolution is sufficiently slow, 
# the system follows the instantaneous ground state and ends up in the 
# solution to the MIS problem.  

# We choose $\Delta_{max} / \Omega_{max} = 2.5$, with $\Omega_{max} = 4 \times 2 \pi$ MHz

# ## Build Pulse Sequence
# We build a discrete adiabatic sweep to analyze the adiabadicity of the two systems by 
# considering and visualizing the minimum gap through their time-dependent energy 
# spectrums.  Because we are considering the weighted MIS problem, we implement individual
# atom detuning with $\Delta(t)_i = w_i \times \Delta(t)$.  We can simulate individual 
# pulse shapes on the emulator.

# return Hamiltonian of system given initial parameters
# return waveforms (for visualization)
function build_adiabatic_sweep(graph, Ω_max::Float64, Δ_max::Float64, t_max::Float64, weights)
    Ω = Waveform(t->Ω_max + (- Ω_max)*t/t_max, duration=t_max)
    Δ = map(1:nv(graph)) do idx
        Waveform(t->weights[idx] * Δ_max * t/t_max, duration=t_max)
    end

    h = XTerm(nv(graph), Ω) - NTerm(nv(graph), Δ)
    return h, Ω, Δ
end

# # Compute MIS Probability and Adiabatic Timescale
# Currently, the emulator does not support directly computing the weighted MIS probability,
# so we use the GraphTensorNetworks package to compute the weighted MIS and corresponding 
# configurations for the two systems.

# We first compute the weighted MIS solution of both original and mapped graphs. The MIS is shown in red.
# Notice that the MIS solution of the mapped graph contains the original MIS (in this case, [1, 3, 4]).

# MIS of the original graph
configs_mapped_0 = solve(IndependentSet(g_0; weights= weight_Δ_0), ConfigsMax())[]
MIS_0 = configs_mapped_0.c[1]
show_graph(g_0; locs = pos_0, vertex_colors=
          [iszero(MIS_0[i]) ? "white" : "red" for i=1:nv(g_0)])
# MIS of mapped graph
configs_mapped_m = solve(IndependentSet(g_m; weights= weight_Δ_M), ConfigsMax())[]
MIS_m = configs_mapped_m.c[1]
show_graph(g_m; locs = simplified_locs, vertex_colors=
          [iszero(MIS_m[i]) ? "white" : "red" for i=1:nv(g_m)])

# In addition to computing the MIS probability, we are also interested in computing 
# the probability that the mapped graph contains the MIS of the original graph 
# We write a function that given a register at a certain time t, and the original graph, 
# we output the MIS probability and the probability the state contains the MIS of the original 
# graph.  These functions assume that all graphs have a unique weighted MIS.  

# this function is a bit cryptic; I'm not sure how to make it easier to understand
function compute_MIS_probability(reg, graph, MIS_original, MIS_mapped)
    v2amp = ThreadsX.map(EaRydCore.ConfigAmplitude(reg)) do (c, amp)
        b = bitarray(c, nv(graph))
        amp_original = 0

        (b[1:length(MIS_original)] == MIS_original) && (amp_original = amp)

        return b, amp, amp_original
    end
    
    return ThreadsX.map(1:2) do k
        ThreadsX.sum(v2amp) do (b, amp, amp_original)
            if k == 1
                if b == MIS_mapped
                    abs2(amp)
                else 
                    zero(real(typeof(amp)))
                end 
            else 
                abs2(amp_original)
            end
                
        end 
    end 
end;

# ## MIS Probability: Original Graph using Independent Set Subspace
# We compute the MIS probability of the original graph as a function of time.  We 
# want to extract the adiabatic timescale $T_{LZ}$ from the Landau-Zener fitting: 
# $1 - P_{MIS} = e^{a - T/T_{LZ}}$.  To do this, we find the first instance time $T^*$
# such that P_{MIS}(T) > 0.9, and continue to run evolutions to 2.5T^* to extract $T_L{LZ}$.
Ω_max = 2 * 4 * pi
Δ_max = 2.5 * Ω_max

# Because the input graph is not a UDG, we must work in the independent set subspace
subspace = independent_set_subspace(g_0)
t_list_o = []
P_MIS_list_o = []

t = 0.1 # length of evolution 
T = 0

while (t < T * 2.5) || (T == 0.0)
    h = build_adiabatic_sweep(g_0, Ω_max, Δ_max, t, weight_Δ_0)[1]
    
    # start evolution at diagonalized state 
    energies, GS = eigsolve(SparseMatrixCSC(h(0.0), subspace), 1, :SR);
    r = RydbergReg(GS[1], subspace);
    
    # run ODE evolution
    prob = ODEEvolution(r, t, h; dt=1e-3, adaptive=false)
    emulate!(prob)
    
    # compute MIS probability
    p = compute_MIS_probability(prob.reg, g_0, MIS_0,MIS_0)[1]

    # find first occurrence p > 0.9
    ((p > 0.9) && (T == 0)) && (T = t)
    
    if (p > 0.9)
        push!(t_list_o, t)
        push!(P_MIS_list_o, p)
    end

    t += 0.01
end

# We can compute the adiabatic timescale by fitting a Landau Zener curve to the 
# MIS probability:

y = broadcast(log, 1 .- P_MIS_list_o)
a, b = linear_fit(t_list_o, y)
T_LZ_0 = -1/b

# plotting settings
gr()
plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)

p = plot(t_list_o, broadcast(log, 1 .- P_MIS_list_o),  palette = :Paired)
scatter!(t_list_o , broadcast(log, 1 .- P_MIS_list_o), title = "MIS Probability", bottom_margin=2cm, left_margin=2cm, right_margin=2cm, label = "Original \$K_{2,4}\$")
plot!(t_list_o, a .+ b .* t_list_o)
xlabel!("\$ T \\ \\ (\\mu s) \$")
ylabel!("log \$ (1 - P_{MIS}) \$")
plot(p, size = (700, 300), fmt = :png)


# We can perform the same analysis with the mapped graph:
subspace_m = independent_set_subspace(g_m)
t_list_m = []
P_MIS_list_m = []
P_MIS_list_m_o = []

t = 0.1
T = 0

while (t < T * 2.5) || (T == 0.0)
    h = build_adiabatic_sweep(g_m, Ω_max, Δ_max, t, weight_Δ_M)[1]
    
    energies, GS = eigsolve(SparseMatrixCSC(h(0.0), subspace_m), 1, :SR);
    r = RydbergReg(GS[1], subspace_m);

    prob = ODEEvolution(r, t, h; dt=1e-3, adaptive=false)
    emulate!(prob)

    p1, p2 = compute_MIS_probability(prob.reg, g_m, MIS_0,MIS_M)

    ((p > 0.9) && (T == 0)) && (T = t)
    
    if (p > 0.9)
        push!(t_list_m, t)
        push!(P_MIS_list_m, p1)
        push!(P_MIS_list_m_o, p2)
    end

    t += 0.01
end

y = broadcast(log, 1 .- P_MIS_list_m)
A, B = linear_fit(t_list_m, y)
T_LZ_m = -1/B

# plot
plot1 = scatter(t_list_m, P_MIS_list_m, legend =:bottomright, title = "Mapped \$K_{2,4}\$: MIS Probability")
scatter!(t_list_m, P_MIS_list_m_o, label = "Mapped \$K_{2,4}\$: Probability contains original MIS ")
xlabel!("\$ T \\ \\ (\\mu s) \$")
ylabel!("\$ P_{MIS} \$")

plot2 = scatter(t_list_m , broadcast(log, 1 .- P_MIS_list_m), title = "MIS Probability", bottom_margin=2cm, left_margin=2cm, right_margin=2cm, label = "Mapped\$K_{2,4}\$")
plot!(t_list_m, A .+ B .* t_list_m)
xlabel!("\$ T \\ \\ (\\mu s) \$")
ylabel!("log \$ (1 - P_{MIS}) \$")

plot(plot1, plot2, layout=(2,1), size = (700, 600), fmt = :png)

























