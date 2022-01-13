# Examples

This directory contains the examples of EaRyd package.

## Setup

Each example contains its own environment, one should start
the corresponding example with its own environment via

```julia
julia --project=examples/<example-name>
```

or from the REPL pkg mode

```julia
] activate examples/<example-name>
```

## Format

The examples are written using [Literate](https://github.com/fredrikekre/Literate.jl).

## Developers Notes

### setup example for development

For developers, it is important to make sure examples are using the
packages in current local directories that is under development. To
do this, just run the following in your terminal

```sh
.ci/run dev <path/to/your/example>
```

this will add the local packages into your example environments.

### create a new example

We enforce each example to have their own environment, an automatic
tool is provided for this purpose, you can create an example folder
with corresponding environment using the following

```sh
.ci/run create <name>
```

which will create `examples/<name>`.
