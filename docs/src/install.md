# [Installation](@id install)

You can copy the following line to your Julia REPL
to install the latest stable version of this package:

```julia
pkg> add Bloqade
```
## Using Bloqade Optimized AMIs from AWS Marketplace
The Bloqade team has implemented 2 dedicated AMIs that can be acquired using AWS Marketplace

- Bloqade AMI with Julia
- Bloqade Optimized Deep Learning AMI with CUDA


### Bloqade AMI with Julia (base image)
The Bloqade AMI contains:
- The latest Julia installation using `juliaup add release` which setup the latest Julia release. For more information please refer to [Juliaup](https://github.com/JuliaLang/juliaup) 
- The latest version of Bloqade 

### Bloqade Optimized Deep Learning AMI with CUDA and Julia
In addition to the content of our base image, this image further contains
- NVIDIA CUDA, cuDNN, NCCL, GPU Drivers, Intel MKL-DNN, Docker, NVIDIA-Docker and EFA support
- Block devices
 /dev/sda1=snap-03d72fbeb983a4663:60:true:gp2
 /dev/sdb=ephemeral0
 /dev/sdc=ephemeral1


### Bloqade AMI(coming soon on AWS Marketplace!)
To use the AMI 
- (Locate the AMIs at AWS Marketplace while selecting the image to start your EC2) "we will add the image ids and arns once all sealed with aws"
- Bloqade.jl/
- follow the [tutorials](https://queracomputing.github.io/Bloqade.jl/dev/) 

## Build System Image to Accelerate Start-up Time

Since Bloqade is a large package, its loading time
and time-to-first-simulation can be very long.
You can build system images to save all the compilation
results in a binary to accelerate its loading/compilation
time. This is useful when you have lots of interactive
programming needs with Bloqade.

To build a system image for your environment, please use
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

## Component Packages

- BloqadeCUDA: CUDA.jl patches for CUDA-based accelerators.
- BloqadeExpr: the interface and expression definition.
- BloqadeKrylov: the Krylov-based solver.
- BloqadeLattices: the lattices definition.
- BloqadeMIS: tools for maximum-independent set.
- BloqadeODE: DiffEq wrapper.
- BloqadeQMC: simulation based on quantum Monte Carlo.
- BloqadeSchema: JSON schema support for Bloqade.
- BloqadeWaveforms: waveform definitions.
- YaoSubspaceArrayReg: the subspace array register for subspace simulation.

## Developing Bloqade

When developing Bloqade, one will need to setup a local environment
that contains all the local changes. To work with the Bloqade repo,
first you need to clone this repo

```sh
# clone this repo
git clone https://github.com/QuEraComputing/Bloqade.jl.git Bloqade
# go into the directory
cd Bloqade
# dev the corresponding environment
.ci/run dev
```

How does this work? The `.ci/run dev` command actually calls the `Pkg.develop`
command from Julia's package manager. Because we want to use the local
changes of the package, one will need to `dev` the corresponding package to 
make the changes happen in your current environment, e.g one will need to `dev` 
the `lib/BloqadeExpr` package to apply changes in `BloqadeExpr` module.

We also provide a convenient tool to setup this more automatically by
looking up dependencies in `lib` in one's `Project.toml` file,

```sh
.ci/run dev <path/to/your/environment>
```

will `dev` all the Bloqade dependencies in your environment.

See also [Modifying A Dependency](https://pkgdocs.julialang.org/v1/getting-started/#Modifying-A-Dependency)
for more detailed explainations.
