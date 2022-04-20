# [Installation](@id install)

You can copy the following line to your Julia REPL
to install the latest stable version of this package.

```julia
pkg> add Bloqade
```

## Build System Image to Accelerate Start-up Time

Since Bloqade is a large package, its loading time
and time-to-first-simulation can be very long.
You can build system images to save all the compilation
results in a binary to accelerate its loading/compilation
time.

### Build System Image Via VSCode Julia Plugin (Recommended)

### Build System Image Via PackageCompiler

## Try the Latest Version of Bloqade

Some users may want to try the latest version of Bloqade for bug fixes, new features, etc. One can use `git` to clone the
repo to try the latest version of the entire package. This
requires one to setup the local project environment via `dev`.
Please refer to the page [Contributing to Bloqade](@ref) for more information.

If you only want to try the latest version of a specific
Bloqade package, just add `#master` after the package name, e.g.

```julia
pkg> add BloqadeExpr#master
```
