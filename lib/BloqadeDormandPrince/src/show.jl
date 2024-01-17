# Forward pretty printing from DormandPrince.jl

function Base.show(io::IO, mime::MIME"text/plain", solver::BloqadeDPSolver)
    printstyled(io, "BloqadeSolver Object", color=:underline)
    println(io, " is using: ")

    show(io, mime, solver.dp_solver)
end