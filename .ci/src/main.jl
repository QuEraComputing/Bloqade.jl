using Pkg
using TOML
using Comonicon
using CoverageTools

root_dir(xs...) = joinpath(dirname(dirname(@__DIR__)), xs...)
function collect_lib(;include_main::Bool=false, excluded_libs=["BloqadePlots"])
    pkgs = Pkg.PackageSpec[]
    for pkg in readdir(root_dir("lib"))
        pkg in excluded_libs && continue
        push!(pkgs, Pkg.PackageSpec(path = root_dir("lib", pkg)))
    end
    include_main && push!(pkgs, Pkg.PackageSpec(path = root_dir()))
    return pkgs
end

"""
create an example.

# Args

- `name`: name of the example.

# Flags

- `-f,--force`: overwrite existing path.
- `--plot`: use `BloqadePlots`.
"""
@cast function create(name::String; force::Bool=false, plot::Bool=false)
    @warn("`.ci/run create` is deprecated, use `.ci/run example create` instead")
    Example.create(name; force, plot)
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
        paths = filter(x->x!=joinpath("lib", "BloqadeCUDA"), paths)
    end
    return test(paths...; coverage)
end

"""
develop the Bloqade components into environment.

# Args

- `path`: path to the environment, relative to the main environment,
    default is the main Bloqade environment.
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
        libs = collect_lib(;include_main=true, excluded_libs=["BloqadeCUDA"])
        Pkg.develop(libs)
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
    libs = readdir(root_dir("lib"))
    d = TOML.parsefile(root_dir(path, "Project.toml"))
    names = [name for name in keys(d["deps"]) if name in libs || name == "Bloqade"]
    paths = map(names) do name
        name == "Bloqade" && return "."
        return root_dir("lib", name)
    end
    pkgs = map(paths) do path
        Pkg.PackageSpec(;path=root_dir(path))
    end
    return pkgs
end

