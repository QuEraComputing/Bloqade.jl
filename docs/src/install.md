# Installation

You can copy the following line to your Julia REPL
to install latest stable version this package.

```julia
pkg> add EaRyd
```

## Build System Image to Accelerate Start-up Time

Since EaRyd is a very large package, its loading time
and time-to-first-emulation can be very long.
You can build system images to save all the compilation
results in a binary to accelerate its loading/compilation
time.

### Build System Image Via VSCode Julia Plugin (Recommended)

### Build System Image Via PackageCompiler

## Try Latest Version of EaRyd

Some users may want to try the latest version of EaRyd for bugfixes, new features etc. One can use `git` to clone the
repo to try the latest version of the entire package. This
requires one to setup the local project environment via `dev`.
Please refer to [#Contributing](@ref) for more reference.

If you only want to try the latest version of a specific
EaRyd package, just add `#master` behind the package name, e.g

```julia
pkg> add EaRydExpr#master
```
