@option struct DiscreteOptions
    progress::Bool = false
    progress_step::Int = 1
    progress_name::String = "emulating"
    normalize_step::Int = 5
    normalize_finally::Bool = true
end

struct KrylovEvolution{S, T <: Real, H <: Hamiltonian}
    reg::S
    durations::Vector{T}
    hamiltonian::H
    options::DiscreteOptions
end

function emulate_step!(prob::KrylovEvolution, t::Real)
    prob.hamiltonian(t)
end
