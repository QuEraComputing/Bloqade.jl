<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/logo.png">
  <img alt="Bloqade Logo">
</picture>
</div>

---

[![CI](https://github.com/QuEraComputing/Bloqade.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/QuEraComputing/Bloqade.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/QuEraComputing/Bloqade.jl/branch/master/graph/badge.svg?token=DYm2XwiTaR)](https://codecov.io/gh/QuEraComputing/Bloqade.jl)

Bloqade is a package developed for quantum computation and quantum simulation based on the neutral-atom architecture with the ability to submit tasks to [QuEra's *Aquila* quantum processor](https://www.quera.com/aquila). Please refer to the [documentation](https://queracomputing.github.io/Bloqade.jl/dev/) page to learn more about Bloqade.

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

## 🚀 Bloqade SDK for Python 🚀

While you're checking out this version of Bloqade, we encourage you to also give [Bloqade for Python](https://github.com/QuEraComputing/bloqade-python) a spin, currently in the Alpha phase of development. Check the documentation to [get started](https://bloqade.quera.com/0.8.0/#how-do-i-get-started). We are in active development of this new package so we appreciate [feedback](https://bloqade.quera.com/0.8.0/contributing/providing-feedback/) during the Alpha phase.

### What About Bloqade.jl?

Bloqade.jl will not go away any time soon because of Bloqade for Python's existence. In fact, plans are underway to allow Bloqade.jl to become a high-performance emulator that the Python version of Bloqade can interface to, giving you unparalleled emulation power to test your ideas on top of a syntax designed to make designing programs for neutral atom hardware easier than ever.

## Citation

If you use Bloqade for a publication, we would kindly ask you to cite our work using the following bibtex citation:

```bibtex
@misc{bloqade2023quera,
  url = {https://github.com/QuEraComputing/Bloqade.jl/},
  title = {Bloqade.jl: {P}ackage for the quantum computation and quantum simulation based on the neutral-atom architecture.},
  year = {2023}
}
```

## License

Apache License 2.0
