using Pkg
using TOML
using Comonicon
using CoverageTools

root_dir(xs...) = joinpath(dirname(dirname(@__DIR__)), xs...)
function collect_lib(; include_main::Bool = false, excluded_libs = [])
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
"""
@cast function create(name::String; force::Bool = false)
    @warn("`.ci/run create` is deprecated, use `.ci/run example create` instead")
    return Example.create(name; force)
end

"""
run tests (in parallel process).

# Args

- `paths`: paths of the packages to run test.

# Flags

- `--coverage`: enable code coverage.
"""
@cast function test(paths::String...; coverage::Bool = false)
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
        return joinpath(path, "src")
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
@cast function testall(; coverage::Bool = false, cpu::Bool = false)
    paths = map(readdir(root_dir("lib"))) do pkg
        return joinpath("lib", pkg)
    end
    push!(paths, ".")
    if cpu
        paths = filter(x -> x != joinpath("lib", "BloqadeCUDA"), paths)
    end
    return test(paths...; coverage)
end

"""
checks for docstrings for packages

# Args

- `paths`: paths of the packages to run test.

"""
@cast function docstring(paths::String...;)
    isempty(paths) && error("expect paths to the package")
    try
        @sync for path in paths
            isfile(joinpath(path, "Manifest.toml")) || dev(path)
            test_cmd = "using Pkg; Pkg.test(;test_args=[\"docstring\"])"
            Threads.@spawn run(`$(Base.julia_exename()) --project=$path -e "$test_cmd"`)
        end
    catch e
        if !(e isa InterruptException)
            rethrow(e)
        end
    end

end


"""
run all the tests for doc strings.

"""
@cast function docstrings()
    paths = map(readdir(root_dir("lib"))) do pkg
        return joinpath("lib", pkg)
    end
    push!(paths, ".")
    return docstring(paths...)
end

"""
develop the Bloqade components into environment.

# Args

- `path`: path to the environment, relative to the main environment,
    default is the main Bloqade environment.

# Flags

- `-v,--verbose`: print the environment status.
"""
@cast function dev(path::String = "."; verbose::Bool = false)
    path = relpath(path, root_dir())

    if path == "."
        Pkg.activate(root_dir(); io = devnull)
        Pkg.develop(collect_lib_deps(path); io = devnull)
    elseif startswith(path, "examples")
        Pkg.activate(root_dir(path); io = devnull)
        Pkg.develop(collect_lib_deps(path); io = devnull)
    elseif startswith(path, "docs") # need all lib packages included
        Pkg.activate(root_dir(path); io = devnull)
        libs = collect_lib(; include_main = true, excluded_libs = ["BloqadeCUDA"])
        Pkg.develop(libs; io = devnull)
    elseif startswith(path, "lib")
        pkgs = collect_lib_deps(path)
        isempty(pkgs) && return
        Pkg.activate(root_dir(path); io = devnull)
        Pkg.develop(pkgs; io = devnull)
    else
        error("invalid path: $(root_dir(path))")
    end

    if verbose
        Pkg.status()
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
        return Pkg.PackageSpec(; path = root_dir(path))
    end
    return pkgs
end
