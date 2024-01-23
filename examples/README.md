# Examples

This directory contains a number of examples of the Bloqade package.

## Generate Jupyter notebooks

you can generate the corresponding jupyter notebook from these examples, e.g

```sh
.ci/run example build 2.adiabatic
```

will build the jupyter notebook at `build` folder. Then you can open and run
this jupyter notebook using the `IJulia` package

```julia
pkg> add IJulia
julia> using IJulia; IJulia.notebook()
```

## Setup

Each example contains its own environment, and one should start
the corresponding example with its own environment via:

```julia
julia --project=examples/<example-name>
```

or from the REPL pkg mode:

```julia
] activate examples/<example-name>
```

## Format

The examples are written using [Literate](https://github.com/fredrikekre/Literate.jl).

## Developers Notes

### Setup example for development

For developers, it is important to make sure examples are using the
packages in current local directories that is under development. To
do this, just run the following in your terminal:

```sh
.ci/run dev <path/to/your/example>
```

This will add the local packages into your example environments. 

If you have any auxiliary files that are accessed from the julia code, put them inside a directory named `data` in your example folder. This directory will be copied to the build directory and can be accessed via a relative path in the script.

### Create a new example

We enforce each example to have their own environment. An automatic
tool is provided for this purpose: you can create an example folder
with the corresponding environment using the following:

```sh
.ci/run create <name>
```

which will create `examples/<name>`. Note that the name should be
unique, and should not contain any whitespaces.
