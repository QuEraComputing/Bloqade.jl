using Pkg
using Literate

root_dir = dirname(@__DIR__)
build_dir = joinpath(root_dir, "build")
example_dir = joinpath(root_dir, "examples")

for each in readdir(example_dir)
    @info "building notebook" path=joinpath(each, "main.jl")
    project_dir = joinpath(example_dir, each)
    isdir(project_dir) || continue
    Pkg.activate(project_dir)
    Literate.notebook(
        joinpath(example_dir, each, "main.jl"),
        project_dir;
        execute=true,
    )
end
