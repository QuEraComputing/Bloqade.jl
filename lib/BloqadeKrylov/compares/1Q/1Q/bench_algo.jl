
using Bloqade
using BloqadeKrylov
using Plots

function H1_SinCos(T, Nt)
    atoms = generate_sites(ChainLattice(), 1, scale = 1)
    wf = Waveform(t->2.2*2π*sin(2π*t), duration = T);
    wf2 = Waveform(t->2.2*2π*cos(2π*t), duration = T);
    h = rydberg_h(atoms; Ω = wf, Δ = wf2)   
    atoms, h
end

function get_exact(T,testing_hamilt)
    NT = 7000
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
    prob = KrylovEvolution(reg, clocks, h)
    emulate!(prob)
    return prob
end

function magnus42(T,NT,testing_hamilt)
    atoms, h = testing_hamilt(T,NT)
    dt = T/NT
    clocks = collect(0:dt:T)
    reg = zero_state(length(atoms))
    prob = CFETEvolution(reg, clocks, h, CFET4_2())
    emulate!(prob)
    return prob
end

function magnus65(T,NT,testing_hamilt)
    atoms, h = testing_hamilt(T,NT)
    dt = T/NT
    clocks = collect(0:dt:T)
    reg = zero_state(length(atoms))
    prob = CFETEvolution(reg, clocks, h, CFET6_5())
    emulate!(prob)
    return prob
end

function magnus811(T,NT,testing_hamilt)
    atoms, h = testing_hamilt(T,NT)
    dt = T/NT
    clocks = collect(0:dt:T)
    reg = zero_state(length(atoms))
    prob = CFETEvolution(reg, clocks, h, CFET8_11())
    emulate!(prob)
    return prob
end

Ns = [2,4,8,16,32,64,128,256,512,1024,2048]
open("Ns","w") do f
    write(f,Ns)
end
Tfull = 1.3
# get exact
res_exact = get_exact(Tfull,H1_SinCos)


## krylov
es = Vector{Float64}(undef,0)
for n in Ns
    res = krylov(Tfull,n,H1_SinCos)
    ## calculate differences:
    e = norm(res_exact.reg.state - res.reg.state)
    push!(es,e)
end
println(es)
open("es_kry","w") do f
    write(f,es)
end

## magnus 42
es = Vector{Float64}(undef,0)
for n in Ns
    res = magnus42(Tfull,n,H1_SinCos)
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
    res = magnus65(Tfull,n,H1_SinCos)
    ## calculate differences:
    e = norm(res_exact.reg.state - res.reg.state)
    push!(es,e)
end
println(es)
open("es_m65","w") do f
    write(f,es)
end

## magnus 811
es = Vector{Float64}(undef,0)
for n in Ns
    res = magnus811(Tfull,n,H1_SinCos)
    ## calculate differences:
    e = norm(res_exact.reg.state - res.reg.state)
    push!(es,e)
end
println(es)
open("es_m811","w") do f
    write(f,es)
end

