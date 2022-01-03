using CSV
using DataFrames
using DelimitedFiles
using EaRydCore
using EaRydODE
using LinearAlgebra
using OrdinaryDiffEq
using Logging: global_logger
using TerminalLoggers: TerminalLogger

include("utils.jl")

function rabi(t, T)
    if t < 0.312
        4/0.312 * t
    elseif 0.312 ≤ t ≤ T + 0.312
        4
    else
        -4/0.312 * (t - T - 0.312 - 0.312)
    end
end

function linear_detuning(t, T)
    p = Polynomials.fit([0.312, 0.312 + T], [15, -11.5], 1)
    if t < 0.312
        15
    elseif 0.312 ≤ t ≤ T + 0.312
        p(t)
    else
        -11.5
    end
end

function rabi_maddie(t, T)
    k = 50
    a = .95
    b = 3.1
    x = t / T

    amplitude = (-1 / (1 + ℯ^(k * (x - a))) ^ b - 1 / (1 + ℯ ^ (-k * (x - (1 - a)))) ^ b + 1) /
                    (-1 / ((1 + ℯ^(k * (1 / 2 - a))) ^ b) - 1 / ((1 + ℯ^(-k * (1 / 2 - (1 - a)))) ^ b) + 1)
    return 4 * amplitude
end

function linear_detuning_maddie(t, T)
    -11 * 2 * (1 / 2 - t / T)
end


function run_pulse(L, R, T_index, graph_index, algo)
    logtimes = LinRange(-2.5, 2, 7)
    times = 2 .^ logtimes
    T = times[T_index]

    lattice = read_lattice(L, graph_index + 1)
    atoms = lattice_to_atoms(lattice)
    graph = unit_disk_graph(atoms, R)
    space = blockade_subspace(graph)
    
    h = RydInteract(atoms) +
        XTerm(length(atoms), t-> 2π * rabi_maddie(t, T)) +
        NTerm(length(atoms), t-> 2π * linear_detuning_maddie(t, T))

    r = zero_state(length(atoms), space)
    emulate!(r, T + 2 * 0.312, h; algo, progress=true, progress_steps=1, reltol=1e-4, abstol=1e-4)

    r.state .= r.state / norm(vec(r.state))
    graph = unit_disk_graph(atoms, 1.5)
    added_rmr = reduced_mean_rydberg(r, graph; add_vertices=true)
    rmr = reduced_mean_rydberg(r, graph; add_vertices=false)
    mis = EaRydCore.exact_solve_mis(graph)
    added_ratio = 1 - added_rmr/mis
    ratio = 1 - rmr/mis
    # save_state(r, L, graph_index, T_index, R)
    # save_ratio(ratio, added_ratio, L, graph_index, T_index, R)
    return ratio
end

run_pulse(5, 1.0, 1, 397)
@time run_pulse(5, 1.5, 7, 397, Tsit5())
@time run_pulse(5, 1.5, 7, 397, Vern8())
@time run_pulse(5, 1.5, 7, 397, Vern6())

@time run_pulse(7, 1.5, 1, 397, Vern6())

logtimes = LinRange(-2.5, 2, 7)
times = 2 .^ logtimes
T = times[1]

lattice = read_lattice(5, 397 + 1)
atoms = lattice_to_atoms(lattice)
graph = unit_disk_graph(atoms, 1.5)
space = blockade_subspace(graph)

h = RydInteract(atoms) +
    XTerm(length(atoms), t-> 2π * rabi_maddie(t, T)) +
    NTerm(length(atoms), t-> 2π * linear_detuning_maddie(t, T))



using SparseArrays
H = SparseMatrixCSC{ComplexF64}(h(1e-3), space)
r = zero_state(length(atoms), space)

eq = EaRydODE.shordinger(h, space)
state = vec(r.state)
dstate = similar(state)

using BenchmarkTools
@benchmark eq($dstate, $state, $h, 2e-3)

@benchmark H * vec(r.state)
@benchmark H * r.state

using CUDA
H = SparseMatrixCSC{ComplexF32}(H)
S = ComplexF32.(vec(r.state))
cuH = CUDA.CUSPARSE.CuSparseMatrixCSC(H)
cuS = cu(S)

C = similar(cuS)
CUDA.allowscalar(false)
typeof(cuH) <: CUDA.CUSPARSE.CuSparseMatrix
typeof(cuS) <: CUDA.CuVector
typeof(C) <: CUDA.CuVector

LinearAlgebra.mul!(C, cuH, cuS)
mul!(C, cuH, cuS)
using BenchmarkTools
@benchmark mul!($C, $cuH, $cuS)

@benchmark CUDA.@sync $cuH * $cuS
@benchmark $H * $S

cuH

using CUDA

r = zero_state(length(atoms), space)
emulate!(r, T + 2 * 0.312, h; algo=Vern8(), progress=true, progress_steps=1)


cuH

using BenchmarkTools

function shordinger(h, s::Subspace; cache=SparseMatrixCSC{ComplexF32}(h(1e-2), s))
    return function equation(dstate, state, h, t)
        H = cu(update_term!(cache, h(t), s))
        H_mul_st = H * state
        @. dstate = -im * H_mul_st
    end
end


CUDA.allowscalar(false)


st = cu(CUDA.rand(ComplexF32, length(r.state)))
prob = ODEProblem(shordinger(h, r.subspace), st, (0.0f0, 0.001f0), h; save_everystep=false, save_start=false, alias_u0=true)
solve(prob, Vern6())
