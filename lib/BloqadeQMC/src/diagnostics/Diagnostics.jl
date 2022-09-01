include("runstats.jl")
include("transition_matrix.jl")

struct Diagnostics{R <: AbstractRunStats, M <: AbstractTransitionMatrix}
    runstats::R
    tmatrix::M

    Diagnostics(runstats::R=NoStats(), tmatrix::M=NoTransitionMatrix()) where {R, M} = new{R, M}(runstats, tmatrix)
end
