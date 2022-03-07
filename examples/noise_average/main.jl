# write your EaRyd example with Literate.jl here
using EaRyd
using CairoMakie
using EaRydPlots
using Random, Distributions



nx, ny = 3, 3
nsites = nx*ny
a= 6.7  # this is the nearest neighbour atom distance 

Random.seed!(123)

σ=0.1  # 0.1 is the atom position flutuation 
# build the atom coordinates 

local d1, d2 = rand(Normal(0, σ), 2)  
atom_coordinate = [(d1, d2)]


for ii = 2: nx
    local d1, d2 = rand(Normal(0, σ ), 2)
    push!(atom_coordinate, ((ii-1)*a+ d1, 0* a+ d2 ))
end

for jj =2: ny
for ii = 1: nx
    local d1, d2 = rand(Normal(0, σ), 2)
    push!(atom_coordinate, ((ii-1)*a +d1, a* (jj-1)+ d2))
end
end

# build the waveform 

# clean part and noise part
total_time = 2.9
Ω_max = 2π * 4.3
U = 2π * 15.0

Ω = Waveform[]
Δ = Waveform[]


for ii=1: nsites

    clock_r = [t for t in 0:1e-2:total_time]
    values_r_o = rand(Normal(0, 0.1), length(clock_r))       
    Ω_f =  piecewise_linear(clocks=clock_r, values= values_r_o)
    Ω_c = piecewise_linear(clocks=[0.0, 0.3, 2.6, total_time], values=[0.0, Ω_max , Ω_max , 0])
    push!(Ω,  Ω_c + Ω_f )
    
    values_r_d = rand(Normal(0, 0.2), length(clock_r))
    Δ_f =  piecewise_linear(clocks=clock_r, values= values_r_d)
    Δ_c = piecewise_linear(clocks=[0.0, 0.3, 2.6, total_time], values=[-U, -U, U , U])

    push!(Δ, Δ_f + Δ_c)
end

h = rydberg_h(atoms; Δ, Ω)


reg = zero_state(9)
prob = ODEEvolution(reg, total_time, h; dt=1e-3, adaptive=false)

densities = []
for info in prob
    push!(densities, [expect(put(nsites, i=>Op.n), info.reg) for i in 1:nsites])
end
D = hcat(densities...)

clocks = [t for t in 0:1e-3:total_time]
heatmap(clocks, 1:9, D'; axis=(xlabel="iterations", ylabel="rydberg density per site"))
