"""
example commands

# Intro

Commands for creating & building examples in examples directory.
"""
@cast module Example

using Pkg
using Comonicon
using ..BloqadeCI: root_dir, collect_lib, dev

function foreach_example(f)
    example_dir = root_dir("examples")
    for name in readdir(example_dir)
        path = joinpath(example_dir, name)
        isdir(path) || continue
        f(path)
    end
    return
end

"""
create an example.

# Args

- `name`: name of the example.

# Flags

- `-f,--force`: overwrite existing path.
"""
@cast function create(name::String; force::Bool = false)
    example_dir = root_dir("examples", name)
    if !force && ispath(example_dir)
        error("$example_dir already exists")
    end
    rm(example_dir; force = true, recursive = true)
    mkpath(example_dir)
    Pkg.activate(example_dir)
    excluded_libs = []
    pkgs = collect_lib(; include_main = true, excluded_libs)
    Pkg.develop(pkgs)
    write(
        joinpath(example_dir, "main.jl"),
        """
# write your Bloqade example with Literate.jl here
""",
    )
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
@cast function build(name::String; build_dir::String = "build", target::String = "notebook", eval::Bool = false)
    ci_dir = root_dir(".ci")
    example_dir = root_dir("examples", name)

    ispath(root_dir(build_dir)) || mkpath(root_dir(build_dir))
    input_file = root_dir("examples", name, "main.jl")
    output_dir = root_dir(build_dir, name)

    if eval
        @info "dev example project" example_dir
        redirect_stdio(stdout = devnull, stdin = devnull, stderr = devnull) do
            return dev(example_dir)
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
        cp(example_dir, output_dir; force = true, follow_symlinks = true)
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
@cast function buildall(; build_dir::String = "build", target::String = "markdown", eval::Bool = false)
    ci_dir = root_dir(".ci")
    example_dir = root_dir("examples")
    script = """
    using Pkg
    using CondaPkg
    using Literate
    for name in readdir(\"$example_dir\")
        project_dir = joinpath(\"$example_dir\", name)
        isdir(project_dir) || continue

        Pkg.activate(project_dir)
        Pkg.instantiate()
        CondaPkg.resolve()

        @info "building" project_dir
        Literate.$target(
            joinpath(project_dir, "main.jl"),
            joinpath(\"$build_dir\", name),
            ;execute=$eval
        )
    end
    """

    # dev the examples first
    # then we run the build in one process
    # so that we can share compile results
    foreach_example() do path
        example_dir = splitpath(path)[end]
        for subdir in readdir(path)
            fullpath = joinpath(path, subdir)
            tutorial_path = 
            isdir(fullpath) && subdir == "data" && run(`cp -r $fullpath $(joinpath(build_dir, example_dir))`)
        end
        return dev(path)
    end
    return run(`$(Base.julia_exename()) --project=$ci_dir -e $script`)
end

end
