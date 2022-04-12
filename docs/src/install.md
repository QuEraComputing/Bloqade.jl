# Installation

You can copy the following line to your Julia REPL
to install latest stable version this package.

```julia
pkg> add Bloqade
```

## Build System Image to Accelerate Start-up Time

Since Bloqade is a very large package, its loading time
and time-to-first-emulation can be very long.
You can build system images to save all the compilation
results in a binary to accelerate its loading/compilation
time.

### Build System Image Via VSCode Julia Plugin (Recommended)

### Build System Image Via PackageCompiler

## Try Latest Version of Bloqade

Some users may want to try the latest version of Bloqade for bugfixes, new features etc. One can use `git` to clone the
repo to try the latest version of the entire package. This
requires one to setup the local project environment via `dev`.
Please refer to [#Contributing](@ref) for more reference.

If you only want to try the latest version of a specific
Bloqade package, just add `#master` behind the package name, e.g

```julia
pkg> add BloqadeExpr#master
```
