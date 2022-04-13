root_dir(args...) = joinpath(homedir(), ".julia", "dev", "Bloqade", args...)
function serve_example(example::String; host::String="0.0.0.0", port::Int=8000)
    # setup environment
    ci_dir = root_dir(".ci")
    docs_dir = root_dir("docs", "build")
    input_file = root_dir("examples", "$(example)", "main.jl")
    output_dir = root_dir("docs", "src", "tutorials")
    project_dir = root_dir("examples", example)
    serve_cmd = """
    using Pkg, FileWatching
    Pkg.activate("$ci_dir")
    using LiveServer, Literate
    Pkg.activate("$project_dir")
    Pkg.instantiate()
    using Bloqade;
    @info "Watching ", "$(project_dir)"
    @async while true
        e = watch_file(joinpath("$project_dir", "main.jl"))
        @info "Event = ", e
        try
            Literate.markdown("$input_file", "$output_dir"; name="$example", execute=true)
            @info "Literate succeed!"
        catch e
            @warn "Literate fail!"
        end
    end
    cd("$docs_dir")
    LiveServer.serve(;
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

serve_example(ARGS[1], port=parse(Int, ARGS[2]))
