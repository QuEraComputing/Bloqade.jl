# [Installation](@id install)

You can copy the following line to your Julia REPL
to install the latest stable version of this package:

```julia
pkg> add Bloqade
```

## Build System Image to Accelerate Start-up Time

Since Bloqade is a large package, its loading time
and time-to-first-simulation can be very long.
You can build system images to save all the compilation
results in a binary to accelerate its loading/compilation
time. This is useful when you have lots of interactive
programming needs with Bloqade.

To build system image for your environment, please use
the [PackageCompiler](https://julialang.github.io/PackageCompiler.jl/dev/)
or use the Julia VSCode plugin's [build system image feature](https://www.julia-vscode.org/docs/stable/userguide/compilesysimage/)

## Try the Latest Version of Bloqade

Some users may want to try the latest version of Bloqade
for bug fixes, new features, etc. One can use `git` to clone the
repo to try the latest version of the entire package. This
requires one to setup the local project environment via `dev`.
Please refer to the page [Contributing to Bloqade](@ref) for more information.

If you only want to try the latest version of a specific
Bloqade package, just add `#master` after the package name, e.g.:

```julia
pkg> add BloqadeExpr#master
```

## Conponent Packages

- BloqadeExpr: the interface and expression definition.
- BloqadeLattices: the lattices definition.
- BloqadeKrylov: the Krylov-based solver.
- BloqadeODE: DiffEq wrapper.
- BloqadeWaveforms: waveform definitions.
- YaoSubspaceArrayReg: the subspace array register for subspace simulation.
- BloqadeCUDA: CUDA.jl patches for CUDA-based accelerators.
- BloqadeMIS: tools for maximum-independent set.

## Developing Bloqade

When developing Bloqade, one will need to setup a local environment
that contains all the local changes. This can be done via the `dev`
command in Julia's Pkg mode. See also [Modifying A Dependency](https://pkgdocs.julialang.org/v1/getting-started/#Modifying-A-Dependency).

One will need to `dev` the corresponding package to make the changes
happen in your current environment, e.g one will need to `dev` the
`lib/BloqadeExpr` package to apply changes in `BloqadeExpr` module.

We also provide a convenient tool to setup this more automatically by
looking up dependencies in `lib` in one's `Project.toml` file,

```sh
.ci/run dev <path/to/your/environment>
```

will `dev` all the Bloqade dependencies in your environment.
