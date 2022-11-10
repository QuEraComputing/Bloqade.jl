struct RydbergPulse{D, B} <: AbstractBlock{D}
    rydberg_hamiltonian::Union{RydbergHamiltonian, RydbergHamiltonian3}
    t_start::Float64
    t_end::Float64
    step::Float64
    function RydbergPulse{D, B}(rydberg_hamiltonian, t_start, t_end, step) where {D, B}
        @assert B isa Union{Type{KrylovEvolution}, Type{SchrodingerProblem}}
        @assert t_start <= t_end
        return new{D, B}(rydberg_hamiltonian, t_start, t_end, step)
    end
end

RydbergPulse(rh::RydbergHamiltonian3, t_start::Real, t_end::Real; backend = KrylovEvolution, step = 1e-2) = RydbergPulse{3, backend}(rh, t_start, t_end, step)
RydbergPulse(rh::RydbergHamiltonian, t_start::Real, t_end::Real; backend = KrylovEvolution, step = 1e-2) = RydbergPulse{2, backend}(rh, t_start, t_end, step)
RydbergPulse(rh::RydbergHamiltonian3, t::Real; backend = KrylovEvolution, step = 1e-2) = RydbergPulse{3, backend}(rh, 0.0, t, step)
RydbergPulse(rh::RydbergHamiltonian, t::Real; backend = KrylovEvolution, step = 1e-2) = RydbergPulse{2, backend}(rh, 0.0, t, step)

function YaoAPI.apply!(reg::AbstractRegister{D}, p::RydbergPulse{D, KrylovEvolution}) where D
    ts = p.t_start:p.step:p.t_end
    prob = KrylovEvolution(reg, ts, p.rydberg_hamiltonian)
    emulate!(prob)
    return reg
end
function YaoAPI.apply!(reg::AbstractRegister{D}, p::RydbergPulse{D, SchrodingerProblem}) where D
    ts = (p.t_start, p.t_end)
    prob = SchrodingerProblem(reg, ts, p.rydberg_hamiltonian)
    emulate!(prob)
    return reg
end
YaoAPI.nqudits(p::RydbergPulse{D, B}) where {D, B} = nqudits(p.rydberg_hamiltonian)
YaoAPI.nlevel(::RydbergPulse{D, B}) where {D, B} = D
YaoAPI.subblocks(p::RydbergPulse{D, B}) where {D, B} = subblocks(p.rydberg_hamiltonian)
function YaoAPI.mat(::Type{T}, p::RydbergPulse{D, B}) where {T, D, B}
    if BloqadeExpr.is_time_dependent(p.rydberg_hamiltonian)
        @warn "Computing the matrix of time-dependent RydbergPulse is slow."
        n = nqudits(p)
        d = nlevel(p)
        mat = Array{T, 2}(undef, d^n, d^n)
        for i in 1:d^n
            st = zeros(T, d^n)
            st[i] = 1
            reg = ArrayReg(st; nlevel = d)
            reg |> p
            mat[:, i] = state(reg)
        end
        return mat
    end
    return exp(T(-im*(p.t_end-p.t_start))*Matrix(YaoAPI.mat(T, p.rydberg_hamiltonian)))
end
YaoAPI.occupied_locs(p::RydbergPulse{D, B}) where {D, B} = tuple(collect(1:nqudits(p))...)
function YaoAPI.print_block(io::IO, p::RydbergPulse{D, B}) where {D, B}
    println(io, "$(nlevel(p))-level Rydberg pulse on $(nqudits(p)) atoms:\nHamiltonian:")
    print_block(io, p.rydberg_hamiltonian)
    println(io, "Time: $(p.t_start) to $(p.t_end)\nSimulated by: $B")
end
Base.show(io::IO, p::RydbergPulse{D, B}) where {D, B} = print_block(io, p)
