# [CPU Acceleration] (@id multithreading)

Bloqade has multithreading-support built-in for faster simulation on multi-core CPUs.

## Backends

This is accomplished through separate Sparse-Matrix dense-Vector (SpMV) *backends* that a user can explicitly choose from to further fine-tune performance. These backends target the CSC format Sparse Matrices that are generated for simulation from Hamiltonian problems. Bloqade also uses other Sparse Matrix formats as well to maximize memory efficiency but with any backend that performs multithreading, the same parallelized functions will be used.

There are three backends to select from:

* `BloqadeExpr` - the default that Bloqade starts with on first installation. Even if you start Julia with multiple threads, when Bloqade has this backend it *will not perform any multithreading*.
* `ThreadedSparseCSR` - Converts CSC matrices to CSR format and uses a simple nested for-loop to perform the SpMV multiplication. The code for this was taken from the [`ThreadedSparseCSR`](https://github.com/BacAmorim/ThreadedSparseCSR.jl) package and migrated into Bloqade to use a more up-to-date version of the [`Polyester`](https://github.com/JuliaSIMD/Polyester.jl) multithreading library, hence the name.
* `ParallelMergeCSR` - Takes the conjugate transpose (adjoint) of CSC matrices and uses QuEra's [`ParallelMergeCSR`](https://github.com/QuEraComputing/ParallelMergeCSR.jl) SpMV package.

!!! compat "Limited Compatibility for ParallelMergeCSR"

    ParallelMergeCSR is currently a *Linux-only* package and may be unreliable/fail to run on other operating systems. You will still be able to take advantage of multithreading through the `ThreadedSparseCSR` backend.

`ThreadedSparseCSR` is ideal for smaller system sizes or if you are performing simulations with the full Hilbert space (in which case the Hamiltonian matrix has a rather even distribution of non-zero elements per row).

`ParallelMergeCSR` is ideal for very large system sizes or if you are performing simulations using the Blockade subspace (@ref subspace) where the imbalance in the number of non-zero entries per row may be larger than if a full Hilbert space simulation was performed. ParallelMergeCSR performs some calculations before the actual matrix-vector multiplication occurs to find the ideal distribution of work across multiple threads which means for smaller system sizes more benefit might be obtained through `ThreadedSparseCSR`.

## Using Multithreading

To enabled multithreading you'll need to first import `BloqadeExpr`.

```julia
julia>using BloqadeExpr
```

You can verify which backend you currently have through:

```julia
BloqadeExpr.backend
```

To set a new backend, pass in the name as a string using the `set_backend` function (you can pass in `"BloqadeExpr"`, `"ThreadedSparseCSR"`, or `"ParallelMergeCSR"`):

```julia
BloqadeExpr.set_backend("ThreadedSparseCSR")
```

You will be prompted to restart the Julia session upon which you should also launch Julia with the desired number of threads:

```
julia -t num_threads
```

All subsequent usage of Bloqade in the Julia session will now take advantage of the full available number of threads.