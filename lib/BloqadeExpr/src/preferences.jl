function set_backend(backend::String)
    if !(backend in ("ParallelMergeCSR", "ThreadedSparseCSR", "BloqadeExpr" ))
        throw(ArgumentError("Invalid Backend provided, please select from \"ParallelMergeCSR\", \"ThreadedSparseCSR\", \"BloqadeExpr\""))
    end

    @set_preferences!("backend" => backend)
    @info("New backend set; restart Julia session for this change to take effect!")
end