# [Contributing to Bloqade](@id contrib)

If you are interested in contributing to this package,
please consider going through this guide to help make your
development workflow as smooth as possible.

Contributing to documentation and unit tests are always great ways to get yourself familiar
with the community and workflows.

## CLI Tool

There is a CLI tool in this repository at `.ci/run` that can help
you simplify your workflow substantially. You can run `.ci/run -h` in your
terminal to see the help message. or run `.ci/run <command> -h`
to see the help message of each command.
Below are some common use cases.

## Documentation

## Setup Documentation

If you are editing the documentation, you can use the `serve` command:

```sh
.ci/run doc serve
```

to serve the documentation locally, and it will automatically update
the served webpage while you are editing. 

If you wish to just build the documentation, you can use `build` command,
which will run the build:

```sh
.ci/run doc build
```

## Light-weight Documentation Build

Due to the "ancient" technology used by Documenter,
it cannot render single page while editing, which
causes the `doc serve` command to be very slow
when you are editing. We provide a light-weight build setup
to workaround this by removing all literate examples
from the documentation. You can enable this by:

```sh
.ci/run doc build --light
.ci/run doc serve --light
```

## Setting Up Environments

The Bloqade package itself is a meta-package that simply re-exports
component packages that live in the `lib` directory. Thus, one will need to
`dev` the corresponding component package to make sure they are
using the `master` branch version while developing. You can always
do this manually in Julia's Pkg mode via the `dev` command. For example, in the
`Bloqade` environment (the `Bloqade/Project.toml` file), one will need
to run the following command:

```julia
pkg> dev lib/BloqadeExpr lib/BloqadeKrylov lib/BloqadeLattices lib/BloqadeMIS lib/BloqadeODE lib/BloqadeWaveforms
```

This can be done automatically using the CLI tool via:

```sh
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
for more detailed explanations.

### Create New Examples

Create a new example project called `my_new_examples` in `examples`
and setup the dependencies of `Bloqade`:

```sh
.ci/run example create my_new_example
```

### Build a Single Example

Build a single example at `build/my_example` to jupyter notebook:

```sh
.ci/run example build my_example
```

### Run Unit Tests

If you are developing unit tests or would like to verify that changes made to any of Bloqade's code have not broken existing functionality, you can run unit tests for a specific sub-package like so:
```sh
.ci/run test path_to_sub_package
```
You can also run ALL tests for the package using 
```sh
.ci/run testall
```
