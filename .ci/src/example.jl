"""
example commands

# Intro

Commands for creating & building examples in examples directory.
"""
@cast module Example

using Pkg
using Comonicon
using ..BloqadeCI: root_dir, collect_lib, dev

"""
create an example.

# Args

- `name`: name of the example.

# Flags

- `-f,--force`: overwrite existing path.
- `--plot`: use `BloqadePlots`.
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
    plot || push!(excluded_libs, "BloqadePlots")
    pkgs = collect_lib(;include_main=true, excluded_libs)
    Pkg.develop(pkgs)
    write(joinpath(example_dir, "main.jl"), """
    # write your Bloqade example with Literate.jl here
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
- `--build-dir=<build dir>`: build directory, default is `build`.

# Flags

- `-e,--eval`: evaluate the Julia code.
"""
@cast function build(name::String; build_dir::String="build", target::String="notebook", eval::Bool=false)
    ci_dir = root_dir(".ci")
    example_dir = root_dir("examples", name)

    ispath(root_dir(build_dir)) || mkpath(root_dir(build_dir))
    input_file = root_dir("examples", name, "main.jl")
    output_dir = root_dir(build_dir, name)

    if eval
        @info "dev example project" example_dir
        redirect_stdio(stdout=devnull, stdin=devnull, stderr=devnull) do
            dev(example_dir)
        end
    end

    setup_env = if eval
        """
        using Pkg
        Pkg.activate(\"$example_dir\")
        Pkg.instantiate()
        """
    else
        ""
    end

    if target == "notebook"
        julia_cmd = """
        using Literate;
        $setup_env
        Literate.notebook(\"$input_file\", \"$output_dir\"; execute=$eval)
        """
        cp(example_dir, output_dir; force=true, follow_symlinks=true)
    elseif target == "markdown"
        julia_cmd = """
        using Literate;
        $setup_env
        Literate.markdown(\"$input_file\", \"$output_dir\"; execute=$eval)
        """
    else
        cmd_error("target $target is not supported")
    end
    run(`$(Base.julia_exename()) --project=$ci_dir -e $julia_cmd`)
    return
end

"""
build all the example in parallel.

# Intro

Similar to `build` but build all the example written Literate
in parallel.

# Options

- `--target=<notebook|markdown>`: build target, either `notebook` or `markdown`.
- `--build-dir=<build dir>`: build directory, default is `build`.

# Flags

- `-e,--eval`: evaluate the Julia code.
"""
@cast function buildall(;build_dir::String="build", target::String="markdown", eval::Bool=false)
    @sync for name in readdir(root_dir("examples"))
        isdir(joinpath(root_dir("examples", name))) || continue
        Threads.@spawn build(name; build_dir, target, eval)
    end
end

end