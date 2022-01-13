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
    include_main && push!(pkgs, Pkg.PackageSpec(path = root_dir("lib")))
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
        Pkg.develop(collect_lib())
    elseif startswith(path, "examples") || startswith(path, "docs")
        Pkg.activate(root_dir(path))
        Pkg.develop(collect_lib(;include_main=true))
    elseif startswith(path, "lib")
        d = TOML.parsefile(joinpath(path, "Project.toml"))
        names = [name for name in keys(d["deps"]) if startswith(name, "EaRyd")]
        isempty(names) && return
        paths = map(names) do name
            name == "EaRyd" && return "."
            return root_dir("lib", name)
        end
        pkgs = map(paths) do path
            Pkg.PackageSpec(;path)
        end
        Pkg.activate(root_dir(path))
        Pkg.develop(pkgs)
    else
        error("invalid path: $(root_dir(path))")
    end
    return
end

"""
create an example.

# Args

- `name`: name of the example.
"""
@cast function create(name::String)
    example_dir = root_dir("examples", name)
    ispath(example_dir) || mkpath(example_dir)
    Pkg.activate(example_dir)
    pkgs = collect_lib()
    push!(pkgs, Pkg.PackageSpec(path = root_dir()))
    Pkg.develop(pkgs)
    write(joinpath(example_dir, "main.jl"), """
    # write EaRyd example with Literate.jl here
    """)
    return
end

"""
documentation commands
"""
@cast module Doc

using Pkg
using Comonicon
using ..EaRydCI: root_dir, dev

"""
build the docs
"""
@cast function build()
    dev("docs")
    Pkg.instantiate()
    docs_dir = root_dir("docs")
    docs_make_jl = root_dir("docs", "make.jl")
    run(`$(Base.julia_exename()) --project=$docs_dir $docs_make_jl`)
end

end

"""
EaRyd CI commands.
"""
@main

end
