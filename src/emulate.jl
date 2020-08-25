function emulate_routine!(r::RydbergReg, ts::Vector, hs::Vector{<: AbstractTerm}, Ks::KrylovSubspace, cache)
    st = vec(r.state)
    for (h, t) in zip(hs, ts)
        update_term!(cache, h, r.subspace)
        # qaoa step
        # NOTE: we share the Krylov subspace here since
        #       the Hamiltonians have the same shape
        arnoldi!(Ks, cache, st)
        expv!(st, -im*t, Ks)
    end
    return r
end

function emulate_routine!(r::Yao.ArrayReg, ts::Vector, hs::Vector{<: AbstractTerm}, Ks::KrylovSubspace, cache)
    st = vec(r.state)
    for (h, t) in zip(hs, ts)
        update_term!(cache, h)
        # qaoa step
        # NOTE: we share the Krylov subspace here since
        #       the Hamiltonians have the same shape
        arnoldi!(Ks, cache, st)
        expv!(st, -im*t, Ks)
    end
    return r
end

"""
    EmulatorCache

Cache type for the emulator.
"""
struct EmulatorCache{K <: KrylovSubspace, C}
    Ks::K
    H::C
end

function EmulatorCache(N::Int, H::AbstractMatrix{Complex{T}}; krylov_niteration=min(30, N)) where T
    Ks = KrylovSubspace{Complex{T}, T}(N, krylov_niteration)
    return EmulatorCache(Ks, H)
end

function EmulatorCache(::Type{T}, H::AbstractTerm, n::Int; kwargs...) where {T <: Real}
    return EmulatorCache(1 << n, SparseMatrixCSC{Complex{T}}(H); kwargs...)
end

function EmulatorCache(::Type{T}, H::AbstractTerm, s::Subspace; kwargs...) where {T <: Real}
    return EmulatorCache(length(s), SparseMatrixCSC{Complex{T}}(H, s); kwargs...)
end

EmulatorCache(H::AbstractTerm, xs...; kwargs...) = EmulatorCache(Float64, H, xs...; kwargs...)

"""
    emulate!(r, ts, hs[; cache=EmulatorCache(ts, hs)])

Emulate the time evolution of a sequence of Hamiltonians `hs` of time length `ts` for each Hamiltonian.
An optional argument can be feeded to preallocate the memory.

    emulate!(r, t::Real, h::AbstractTerm[; algo=Tsit5(), kwargs...])

Emulate the contiguous time evolution of given Hamiltonian `h` that has dependency on `t`. The default algorithm
we use is `Tsit5`. Available `kwargs` can be found in [Common Solver Options](https://diffeq.sciml.ai/latest/basics/common_solver_opts/).
Available algorithms can be found at [Full List of Methods](https://diffeq.sciml.ai/latest/solvers/ode_solve/#Full-List-of-Methods).
"""
function emulate! end

function emulate!(r::Yao.ArrayReg, ts::Vector{T}, hs::Vector{<:AbstractTerm}; cache=EmulatorCache(T, first(hs), nsites(first(hs)))) where {T <: Real}
    emulate_routine!(r, ts, hs, cache.Ks, cache.H)
    return r
end

function emulate!(r::RydbergReg, ts::Vector{T}, hs::Vector{<:AbstractTerm}; cache=EmulatorCache(T, first(hs), r.subspace)) where {T <: Real}
    emulate_routine!(r, ts, hs, cache.Ks, cache.H)
    return r
end


# contiguous
"""
    shordinger(s::Subspace)

Create a Shordinger equation in given subspace `s`.
"""
function shordinger(h::AbstractTerm, s::Subspace)
    H = SparseMatrixCSC(h(0.1), s)

    return function equation(dstate, state, h, t)
        dstate .= -im .* update_term!(H, h(t), s) * state
    end
end

"""
    shordinger()

Create a Shordinger equation.
"""
function shordinger(h::AbstractTerm)
    H = SparseMatrixCSC(h(0.1))

    return function equation(dstate, state, h, t)
        dstate .= -im .* update_term!(H, h(t)) * state
    end
end

function emulate!(r::Yao.ArrayReg, t::Real, h::AbstractTerm; algo=Tsit5(), kwargs...)
    prob = ODEProblem(shordinger(h), vec(r.state), (zero(t), t), h; save_everystep=false, save_start=false, alias_u0=true, kwargs...)
    result = solve(prob, algo)
    return r
end

function emulate!(r::RydbergReg, t::Real, h::AbstractTerm; algo=Tsit5(), kwargs...)
    prob = ODEProblem(shordinger(h, r.subspace), vec(r.state), (zero(t), t), h; save_everystep=false, save_start=false, alias_u0=true, kwargs...)
    result = solve(prob, algo)
    return r
end
