struct RydbergPulse{D} <: AbstractBlock{D}
    rydberg_hamiltonian::Union{RydbergHamiltonian, RydbergHamiltonian3}
    t_start::Float64
    t_end::Float64
    backend::Union{Type{KrylovEvolution}, Type{SchrodingerProblem}}
    step::Float64
end

RydbergPulse(rh::RydbergHamiltonian3, t_start::Real, t_end::Real; backend = KrylovEvolution, step = 1e-2) = RydbergPulse{3}(rh, t_start, t_end, backend, step)
RydbergPulse(rh::RydbergHamiltonian, t_start::Real, t_end::Real; backend = KrylovEvolution, step = 1e-2) = RydbergPulse{2}(rh, t_start, t_end, backend, step)
RydbergPulse(rh::RydbergHamiltonian3, t::Real; backend = KrylovEvolution, step = 1e-2) = RydbergPulse{3}(rh, 0.0, t, backend, step)
RydbergPulse(rh::RydbergHamiltonian, t::Real; backend = KrylovEvolution, step = 1e-2) = RydbergPulse{2}(rh, 0.0, t, backend, step)

function YaoAPI.apply!(reg::AbstractRegister{D}, p::RydbergPulse{D}) where D
    evo = p.backend
    if evo === KrylovEvolution
        @assert p.t_start <= p.t_end
        ts = p.t_start:p.step:p.t_end
    else
        ts = (p.t_start, p.t_end)
    end
    prob = evo(reg, ts, p.rydberg_hamiltonian)
    emulate!(prob)
    return reg
end

Base.:(|>)(reg::AbstractRegister{D}, p::RydbergPulse{D}) where D = apply!(reg, p)
YaoAPI.nqudits(p::RydbergPulse{D}) where D = nqudits(p.rydberg_hamiltonian)
YaoAPI.nlevel(::RydbergPulse{D}) where D = D
function YaoAPI.mat(p::RydbergPulse{D}) where D 
    if BloqadeExpr.is_time_dependent(p.rydberg_hamiltonian)
        error("Cannot get the matrix representation of a time-dependent Hamiltonian evolution")
    end
    return exp(-im*(p.t_end-p.t_start)*Matrix(p.rydberg_hamiltonian))
end
YaoAPI.occupied_locs(p::RydbergPulse{D}) where D = 1:nqudits(p)
function YaoAPI.print_block(io::IO, p::RydbergPulse{D}) where D
    print_block(io, p.rydberg_hamiltonian)
    println(io, "time: $(p.t_start) to $(p.t_end)\t simulated by: $(p.backend)")
end
Base.show(io::IO, p::RydbergPulse{D}) where D = print_block(io, p)
