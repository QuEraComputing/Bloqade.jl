
using Bloqade
using BloqadeKrylov
using Plots



function H2_Square(T,Nt)
    nx, ny = 3, 3
    nsites = nx * ny
    atoms = generate_sites(SquareLattice(), nx, ny, scale = 6.7)
    
    Ω_max = 2π * 4.3
    Ω = piecewise_linear(clocks = [0.0/2.9, 0.3/2.9, 2.6/2.9, 1.0].*T, values = [0.0, Ω_max, Ω_max, 0]);

    U = 2π * 15.0
    Δ = piecewise_linear(clocks = [0.0/2.9, 0.3/2.9, 2.6/2.9, 1.0].*T, values = [-U, -U, U, U]);

    h = rydberg_h(atoms; Δ, Ω)
    atoms, h
end

function get_exact(T,testing_hamilt)
    NT = 10000
    atoms, h = testing_hamilt(T,NT)
    dt = T/NT
    clocks = collect(0:dt:T)
    odereg = zero_state(length(atoms))
    ODEprob = SchrodingerProblem(odereg,T,h)
    emulate!(ODEprob)
    return ODEprob
end

function krylov(T,NT,testing_hamilt)
    atoms, h = testing_hamilt(T,NT)
    dt = T/NT
    clocks = collect(0:dt:T)
    reg = zero_state(length(atoms))
    prob = KrylovEvolution(reg, clocks, h, tol=1.0e-12)
    emulate!(prob)
    return prob
end

function magnus42(T,NT,testing_hamilt)
    atoms, h = testing_hamilt(T,NT)
    dt = T/NT
    clocks = collect(0:dt:T)
    reg = zero_state(length(atoms))
    prob = CFETEvolution(reg, clocks, h, CFET4_2(),tol=1.0e-12)
    emulate!(prob)
    return prob
end

function magnus65(T,NT,testing_hamilt)
    atoms, h = testing_hamilt(T,NT)
    dt = T/NT
    clocks = collect(0:dt:T)
    reg = zero_state(length(atoms))
    prob = CFETEvolution(reg, clocks, h, CFET6_5(),tol=1.0e-12)
    emulate!(prob)
    return prob
end

function magnus811(T,NT,testing_hamilt)
    atoms, h = testing_hamilt(T,NT)
    dt = T/NT
    clocks = collect(0:dt:T)
    reg = zero_state(length(atoms))
    prob = CFETEvolution(reg, clocks, h, CFET8_11(),tol=1.0e-12)
    emulate!(prob)
    return prob
end

Ns = [2,4,8,16,32,64,128,256,512,1024,2048,4069,8192]
open("Ns","w") do f
    write(f,Ns)
end
Tfull = 2.9

# get exact
res_exact = get_exact(Tfull,H2_Square)


## magnus 42
es = Vector{Float64}(undef,0)
for n in Ns
    res = magnus42(Tfull,n,H2_Square)
    ## calculate differences:
    e = norm(res_exact.reg.state - res.reg.state)
    push!(es,e)
end
println(es)
open("es_m42","w") do f
    write(f,es)
end

## magnus 65
es = Vector{Float64}(undef,0)
for n in Ns
    res = magnus65(Tfull,n,H2_Square)
    ## calculate differences:
    e = norm(res_exact.reg.state - res.reg.state)
    push!(es,e)
end
println(es)
open("es_m65","w") do f
    write(f,es)
end

## krylov
es = Vector{Float64}(undef,0)
for n in Ns
    res = krylov(Tfull,n,H2_Square)
    ## calculate differences:
    e = norm(res_exact.reg.state - res.reg.state)
    push!(es,e)
end
println(es)
open("es_kry","w") do f
    write(f,es)
end



## magnus 811
es = Vector{Float64}(undef,0)
for n in Ns
    res = magnus811(Tfull,n,H2_Square)
    ## calculate differences:
    e = norm(res_exact.reg.state - res.reg.state)
    push!(es,e)
end
println(es)
open("es_m811","w") do f
    write(f,es)
end

