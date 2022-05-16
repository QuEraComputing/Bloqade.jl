using Test
using Random
using LinearAlgebra
using YaoArrayRegister
using BitBasis

Random.seed!(8)
test_subspace = Subspace(5, sort(randperm(32)[1:32]) .- 1)
# prepare a zero state
r = zero_state(test_subspace)
sample1 = measure!(r)
@test sample1 == zero(BitStr64{5})
# 1. sampling
r = rand_state(test_subspace)
samples = measure(r; nshots = 10000)
@test samples isa Vector{<:BitStr}
