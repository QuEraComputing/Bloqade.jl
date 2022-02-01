module EaRydCI

using Pkg
using TOML
using Comonicon
using CoverageTools

root_dir(xs...) = joinpath(dirname(dirname(@__DIR__)), xs...)
function collect_lib(;include_main::Bool=false, excluded_libs=["EaRydPlots"])
    pkgs = Pkg.PackageSpec[]
    for pkg in readdir(root_dir("lib"))
        pkg in excluded_libs && continue
        push!(pkgs, Pkg.PackageSpec(path = root_dir("lib", pkg)))
    end
    include_main && push!(pkgs, Pkg.PackageSpec(path = root_dir()))
    return pkgs
end

"""
run tests (in parallel process).

# Args

- `paths`: paths of the packages to run test.

# Flags

- `--coverage`: enable code coverage.
"""
@cast function test(paths::String...; coverage::Bool=false)
    isempty(paths) && error("expect paths to the package")
    try
        @sync for path in paths
            isfile(joinpath(path, "Manifest.toml")) || dev(path)
            test_cmd = "using Pkg; Pkg.test(;coverage=$coverage)"
            Threads.@spawn run(`$(Base.julia_exename()) --project=$path -e "$test_cmd"`)
        end
    catch e
        if !(e isa InterruptException)
            rethrow(e)
        end
    end

    coverage || return
    
    dirs = map(paths) do path
        joinpath(path, "src")
    end
    pfs = mapreduce(process_folder, vcat, dirs)
    LCOV.writefile(root_dir("lcov.info"), pfs)
    return
end

"""
run all the tests (in parallel).

# Flags

- `--coverage`: enable code coverage.
- `--cpu`: filter cpu tests.
"""
@cast function testall(;coverage::Bool=false, cpu::Bool=false)
    paths = map(readdir(root_dir("lib"))) do pkg
        joinpath("lib", pkg)
    end
    push!(paths, ".")
    if cpu
        paths = filter(x->x!=joinpath("lib", "EaRydCUDA"), paths)
    end
    return test(paths...; coverage)
end

"""
develop the EaRyd components into environment.

# Args

- `path`: path to the environment, relative to the main environment,
    default is the main EaRyd environment.
"""
@cast function dev(path::String=".")
    path = relpath(path, root_dir())

    if path == "."
        Pkg.activate(root_dir())
        Pkg.develop(collect_lib_deps(path))
    elseif startswith(path, "examples")
        Pkg.activate(root_dir(path))
        Pkg.develop(collect_lib_deps(path))
    elseif startswith(path, "docs") # need all lib packages included
        Pkg.activate(root_dir(path))
        Pkg.develop(collect_lib_deps(path))
    elseif startswith(path, "lib")
        pkgs = collect_lib_deps(path)
        isempty(pkgs) && return
        Pkg.activate(root_dir(path))
        Pkg.develop(pkgs)
    else
        error("invalid path: $(root_dir(path))")
    end
    return
end

function collect_lib_deps(path::String)
    d = TOML.parsefile(root_dir(path, "Project.toml"))
    names = [name for name in keys(d["deps"]) if startswith(name, "EaRyd")]
    paths = map(names) do name
        name == "EaRyd" && return "."
        return root_dir("lib", name)
    end
    pkgs = map(paths) do path
        Pkg.PackageSpec(;path=root_dir(path))
    end
    return pkgs
end

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
documentation commands
"""
@cast module Doc

using Pkg
using Comonicon
using LiveServer
using ..EaRydCI: root_dir, dev

@cast function serve(;host::String="0.0.0.0", port::Int=8000)
    # setup environment
    dev("docs")
    docs_dir = root_dir("docs")
    julia_cmd = "using Pkg; Pkg.instantiate()"
    run(`$(Base.julia_exename()) --project=$docs_dir -e $julia_cmd`)


    docs_dir = root_dir(".ci")
    serve_cmd = """
    using LiveServer;
    LiveServer.servedocs(;
        doc_env=true,
        skip_dirs=[
            joinpath("docs", "src", "assets"),
            joinpath("docs", "src", "tutorials"),
        ],
        literate="examples",
        host=\"$host\",
        port=$port,
    )
    """
    try
        run(`$(Base.julia_exename()) --project=$(root_dir(".ci")) -e $serve_cmd`)
    catch e
        if e isa InterruptException
            return
        else
            rethrow(e)
        end
    end
    return
end

"""
build the docs
"""
@cast function build()
    dev("docs")
    docs_dir = root_dir("docs")
    docs_make_jl = root_dir("docs", "make.jl")
    julia_cmd = "using Pkg; Pkg.instantiate()"
    
    for each in readdir(root_dir("examples"))
        example_dir = root_dir("examples", each)
        isdir(example_dir) || continue
        dev(example_dir)
    end
    run(`$(Base.julia_exename()) --project=$docs_dir -e $julia_cmd`)
    run(`$(Base.julia_exename()) --project=$docs_dir $docs_make_jl`)
end

end

"""
EaRyd CI commands.
"""
@main

end
