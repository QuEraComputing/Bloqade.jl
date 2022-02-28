"""
example commands
"""
@cast module Example

using Pkg
using Comonicon
using ..EaRydCI: root_dir, collect_lib

"""
create an example.

# Args

- `name`: name of the example.

# Flags

- `-f,--force`: overwrite existing path.
- `--plot`: use `EaRydPlots`.
"""
@cast function create(name::String; force::Bool=false, plot::Bool=false)
    example_dir = root_dir("examples", name)
    if !force && ispath(example_dir)
        error("$example_dir already exists")
    end
    rm(example_dir;force=true, recursive=true)
    mkpath(example_dir)
    Pkg.activate(example_dir)
    excluded_libs = []
    plot || push!(excluded_libs, "EaRydPlots")
    pkgs = collect_lib(;include_main=true, excluded_libs)
    Pkg.develop(pkgs)
    write(joinpath(example_dir, "main.jl"), """
    # write your EaRyd example with Literate.jl here
    """)
    return
end

"""
build the example.

# Intro

Build a single example to `build` directory, by default
it generates the jupyter notebook of the corresponding
Literate script and copy all other files to the build directory.

# Args

- `name`: name of the example project you would like to build.

# Options

- `--target=<notebook|markdown>`: build target, either `notebook` or `markdown`.

# Flags

- `-e,--eval`: evaluate the Julia code.
"""
@cast function build(name::String; target::String="notebook", eval::Bool=false)
    ci_dir = root_dir(".ci")
    example_dir = root_dir("examples", name)

    ispath(root_dir("build")) || mkpath(root_dir("build"))
    input_file = root_dir("examples", name, "main.jl")
    output_dir = root_dir("build", name)
    if target == "notebook"
        julia_cmd = """
        using Pkg, Literate;
        Pkg.activate(\"$example_dir\")
        Literate.notebook(\"$input_file\", \"$output_dir\"; execute=$eval)
        """
        cp(example_dir, output_dir; force=true, follow_symlinks=true)
    elseif target == "markdown"
        julia_cmd = """
        using Pkg, Literate;
        Pkg.activate(\"$example_dir\")
        Literate.markdown(\"$input_file\", \"$output_dir\"; execute=$eval)
        """
    else
        cmd_error("target $target is not supported")
    end
    run(`$(Base.julia_exename()) --project=$ci_dir -e $julia_cmd`)
    return
end

end