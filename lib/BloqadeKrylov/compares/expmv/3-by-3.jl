include("tools.jl")

## define problem-1
## 1) real time
prob_name = "3-by-3"

function H2_Square(T,ascale)
    nx, ny = 3, 3
    nsites = nx * ny
    atoms = generate_sites(SquareLattice(), nx, ny, scale = ascale)
    
    Ω_max = 2π * 4.3
    Ω = piecewise_linear(clocks = [0.0/T, 0.3/T, 2.6/T, 1.0].*T, values = [0.0, Ω_max, Ω_max, 0]);

    U = 2π * 15.0
    Δ = piecewise_linear(clocks = [0.0/T, 0.3/T, 2.6/T, 1.0].*T, values = [-U, -U, U, U]);

    h = rydberg_h(atoms; Δ, Ω)
    atoms, h
end

nsites = 9

T = 2.9
Nt = 1000
scale = 3.0

atoms, h = H2_Square(T,scale)
dt = T/Nt
tspan = dt:dt:T

H = BloqadeExpr.Hamiltonian(ComplexF64, h)

# expmv
reg = zero_state(nsites)
nmul_old = testing_expmv!(H, tspan, reg)

# expmv
reg1 = zero_state(nsites)
#tspan = dt:dt:2dt
nmul_overhead, nmul_impl = testing_expm_multiply!(H, tspan, reg1)


#checking
#print(statevec(reg1))
#print(statevec(reg))
println(statevec(reg) .≈ statevec(reg1))



open("$prob_name.$Nt._$(scale)_.exmp","w") do f
    write(f,nmul_old)
end

open("$prob_name.$Nt._$(scale)_.exm_mply.ovh","w") do f
    write(f,nmul_overhead)
end
open("$prob_name.$Nt._$(scale)_.exm_mply.impl","w") do f
    write(f,nmul_impl)
end

open("$prob_name.$Nt._$(scale)_.ts","w") do f
    write(f,collect(tspan))
end