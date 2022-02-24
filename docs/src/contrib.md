# Contributing

If you are interested in contributing to this package,
please consider going through this guide to help you smooth your
developing workflow.

Contributing to documentation is always a good start to get familiar
with the community and workflows.

## Documentation

## Setting Up Environments

## Components

## The CLI Tool

There is a CLI tool in this repository at `.ci/run` that can help
you simplify the workflow a lot. You can run `.ci/run -h` in your
terminal to see the help message. or run `.ci/run <command> -h`
to see the help message of each command.
Here are some common examples of it.

### Create New Examples

create a new example project called `my_new_examples` in `examples`
and setup the dependencies of `EaRyd`.

```sh
.ci/run example create my_new_example
```

### Build Single Example

build a single example at `build/my_example` to jupyter notebook.

```sh
.ci/run example build my_example
```

