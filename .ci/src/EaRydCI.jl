module EaRydCI

using Pkg
using TOML
using Comonicon

root_dir(xs...) = joinpath(dirname(dirname(@__DIR__)), xs...)
function collect_lib()
    return map(readdir(root_dir("lib"))) do pkg
        Pkg.PackageSpec(path = root_dir("lib", pkg))
    end
end

@cast function test(paths::String...; coverage::Bool=false)
    isempty(paths) && error("expect paths to the package")
    paths = collect(paths)
    names = map(paths) do pkg
        d = TOML.parsefile(joinpath(pkg, "Project.toml"))
        d["name"]
    end

    isfile(root_dir("Manifest.toml")) || dev(".") # ensure all lib is developed
    test_cmd = "using Pkg; Pkg.test([$(join(repr.(names), ", "))];coverage=$coverage)"
    run(`$(Base.julia_exename()) --project=$(root_dir()) -e "$test_cmd"`)
    return
end

"""
develop the EaRyd components into environment.

# Args

- `path`: path to the environment, relative to the main environment,
    default is the main EaRyd environment.
"""
@cast function dev(path::String=".")
    pkgs = collect_lib()
    path = relpath(path, root_dir())

    if path == "."
        Pkg.activate(root_dir())
        Pkg.develop(pkgs)
    elseif startswith(path, "examples") || startswith(path, "docs")
        Pkg.activate(root_dir(path))
        push!(pkgs, Pkg.PackageSpec(path = root_dir()))
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
