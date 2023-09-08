include("tools.jl")

## define problem-1
## 1) real time
prob_name = "chain_sin_cos"


nsites = 10
dt = 1e-3
tspan = 0:dt:1.0
atoms = generate_sites(ChainLattice(), nsites, scale = 6)
h = rydberg_h(atoms; Ω = sin, ϕ = cos)

H = BloqadeExpr.Hamiltonian(ComplexF64, h)

# expmv
reg = zero_state(nsites)
nmul_old = testing_expmv!(H, tspan, reg)

# expmv
reg1 = zero_state(nsites)
nmul_overhead, nmul_impl = testing_expm_multiply!(H, tspan, reg1)

#checking
println(statevec(reg) ≈ statevec(reg1))


open("$prob_name.$nsite","w") do f
    write(f,es)
end

