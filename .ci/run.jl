using Pkg

function main()
    root_directory = dirname(@__DIR__)
    package_names = readdir(joinpath(root_directory, "lib"))

    help = """
    EaRyd CI Test Manager

    Usage

        run.jl <command>

    Commands

    dev                 develop the packages in current environment
    test <device>       test the packages on cpu or cuda.
    doc                 
    """

    length(ARGS) > 0 || return print(help)
    if "dev" == ARGS[1]
        packages = map(package_names) do pkg
            Pkg.PackageSpec(path = joinpath(root_directory, "lib", pkg))
        end
        Pkg.develop(packages)
    elseif "test" == ARGS[1]
        length(ARGS) == 2 || return print(help)
        device = ARGS[2]
        if device == "cuda"
        elseif device == "cpu"
            package_names = filter(!startswith("Cu"), package_names)
        else
            error("invalid device")
        end
        Pkg.test([package_names..., "EaRyd"]; coverage=true)
    elseif "doc" == ARGS[1]
        packages = map(package_names) do pkg
            Pkg.PackageSpec(path = joinpath(root_directory, "lib", pkg))
        end
        push!(packages, Pkg.PackageSpec(path = root_directory))
        Pkg.develop(packages)
        Pkg.instantiate()
    end
end

main()
