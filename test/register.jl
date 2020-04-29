using Test
using RydbergEmulator
using Yao

subspace = [0, 1, 4, 8]
raw_st = zeros(ComplexF64, length(subspace), 5)
raw_st[1, :] .= 1
@test state(RydbergEmulator.zero_state(10, subspace; nbatch=5)) ≈ raw_st
