using Pkg

function main()
    root_directory = dirname(@__DIR__)
    @info root_directory
    package_names = readdir(joinpath(root_directory, "lib"))
    @info "package names" package_names

    help = """
    EaRyd CI Test Manager

    Usage

        run.jl <command>

    Commands

    dev                         develop the packages in current environment
    test <device> [--cuda]      test the packages on cpu or cuda.
    doc                         dev all doc dependencies
    example                     dev all example dependencies at current environment
    """

    length(ARGS) > 0 || return print(help)
    if "dev" == ARGS[1]
        packages = map(package_names) do pkg
            Pkg.PackageSpec(path = joinpath(root_directory, "lib", pkg))
        end
        Pkg.develop(packages)
    elseif "test" == ARGS[1]
        length(ARGS) ≤ 3 || return print(help)

        if length(ARGS) == 2 && "--cpu" in ARGS # test all on cpu
            package_names = filter(!endswith("CUDA"), package_names)
        elseif length(ARGS) == 2 && "--cuda" in ARGS # test all on cuda
            push!(package_names, "EaRyd")
        elseif length(ARGS) == 2
            ARGS[2] == "EaRyd" || ARGS[2] in package_names || return print(help)
            package_names = [ARGS[2]]
        else
            return print(help)
        end
        Pkg.test(package_names; coverage=true)
    elseif "doc" == ARGS[1]
        packages = map(package_names) do pkg
            Pkg.PackageSpec(path = joinpath(root_directory, "lib", pkg))
        end
        push!(packages, Pkg.PackageSpec(path = root_directory))
        Pkg.develop(packages)
        Pkg.instantiate()
    elseif "example" == ARGS[1]
        packages = map(package_names) do pkg
            Pkg.PackageSpec(path = joinpath(root_directory, "lib", pkg))
        end
        push!(packages, Pkg.PackageSpec(path = root_directory))
        Pkg.develop(packages)
    end
end

main()
