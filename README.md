<div align="center"> <img
src="docs/src/assets/logo.png"
alt="Bloqade Logo" width="400"></img>
</div>

---

[![Build Status](https://github.com/QuEraComputing/Bloqade.jl/workflows/CI/badge.svg)](https://github.com/QuEraComputing/Bloqade.jl/actions)
[![codecov](https://codecov.io/gh/QuEraComputing/Bloqade.jl/branch/master/graph/badge.svg?token=DYm2XwiTaR)](https://codecov.io/gh/QuEraComputing/Bloqade.jl)

**If you have the following issue installing Bloqade**

```
  ✗ Bloqade
  243 dependencies successfully precompiled in 550 seconds (16 already precompiled)
  1 dependency errored. To see a full report either run `import Pkg; Pkg.precompile()` or load the package
     Testing Running tests...
    CondaPkg Found dependencies: /home/runner/work/Bloqade.jl/Bloqade.jl/CondaPkg.toml
    CondaPkg Found dependencies: /home/runner/.julia/packages/PythonCall/Td3SH/CondaPkg.toml
    CondaPkg Resolving changes
             + libstdcxx-ng
             + matplotlib
             + python
 Downloading artifact: micromamba-0.27.0
 Downloading artifact: micromamba-0.27.0
┌ Error: Hash Mismatch!
│   Expected sha256:   4adbf3091a4159af2c48264a8e32ecb98147b0e3f200601f384f8f53a6910ca2
│   Calculated sha256: 1bb0c8896927a64a6d73a33fa08a915c22c57b240db92e2d6595b6741f509ed0
└ @ Pkg.PlatformEngines /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6/Pkg/src/PlatformEngines.jl:629
ERROR: LoadError: InitError: Unable to automatically install 'micromamba-0.27.0' from '/home/runner/.julia/packages/MicroMamba/rCGZ4/Artifacts.toml'
Stacktrace:
  [1] error(s::String)
```

**please consider run the following command instead**

```julia
pkg> add https://github.com/kshyatt/MicroMamba.jl.git#ksh/artifactfix

pkg> add Bloqade
```

Bloqade is a package developed for quantum computation and quantum simulation based on the neutral-atom architecture. Please refer to the [documentation](https://queracomputing.github.io/Bloqade.jl/dev/) page to learn more about Bloqade.

**Bloqade is currently under public release beta. High-level APIs are stable, but please expect some further exploration and rough edges.**

## Installation

<p>
Bloqade is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package. To install Bloqade,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode, and then type the following command:
</p>

For stable release:

```julia
pkg> add Bloqade
```

For current master:

```julia
pkg> add Bloqade#master
```

## Community 

You can join QuEra's Slack workspace with this [link](https://join.slack.com/t/querapublic/shared_invite/zt-1d5jjy2kl-_BxvXJQ4_xs6ZoUclQOTJg). Join the `#bloqade` channel to discuss anything related to Bloqade.

## License

Apache License 2.0
