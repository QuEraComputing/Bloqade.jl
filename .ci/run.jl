using Pkg
root_directory = dirname(@__DIR__)

package_names = readdir(joinpath(root_directory, "lib"))
device = get(ENV, "EARYD_TEST_DEVICE", "cpu")
if device == "cuda"
elseif device == "cpu"
    package_names = filter(!startswith("Cu"), package_names)
end

packages = map(package_names) do pkg
    Pkg.PackageSpec(path = joinpath(root_directory, "lib", pkg))
end

help = """
EaRyd CI Test Manager

Usage

    run.jl <command>

Commands

dev         develop the packages in current environment
test        test the packages
"""

length(ARGS) == 1 || print(help)

if "dev" == only(ARGS)
    Pkg.develop(packages)
elseif "test" == only(ARGS)
    Pkg.test([package_names..., "EaRyd"]; coverage=true)
end
