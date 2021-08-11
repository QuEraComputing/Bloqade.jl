# Copyright 2020 QuEra Computing Inc. All rights reserved.

module ContinuousEmulator

using Reexport
@reexport using RydbergEmulator
using OrdinaryDiffEq
using Yao
using SparseArrays
using RydbergEmulator: AbstractTerm

"""
    shordinger(h::AbstractTerm, s::Subspace; cache=SparseMatrixCSC(h(1e-2), s))

Create a Shordinger equation in given subspace `s`.
"""
function shordinger(h::AbstractTerm, s::Subspace; cache=SparseMatrixCSC(h(1e-2), s))
    return function equation(dstate, state, h, t)
        dstate .= -im .* update_term!(cache, h(t), s) * state
    end
end

"""
    shordinger(h::AbstractTerm; cache=SparseMatrixCSC(h(1e-2)))

Create a Shordinger equation.
"""
function shordinger(h::AbstractTerm; cache=SparseMatrixCSC(h(1e-2)))
    return function equation(dstate, state, h, t)
        dstate .= -im .* update_term!(cache, h(t)) * state
    end
end

function RydbergEmulator.emulate!(r::Yao.ArrayReg, t::Real, h::AbstractTerm; algo=lsoda(), cache=SparseMatrixCSC(h(1e-2)), kwargs...)
    prob = ODEProblem(shordinger(h; cache=cache), vec(r.state), (zero(t), t), h; save_everystep=false, save_start=false, alias_u0=true, kwargs...)
    result = solve(prob, algo)
    return r
end

function RydbergEmulator.emulate!(r::RydbergReg, t::Real, h::AbstractTerm; algo=lsoda(), cache=SparseMatrixCSC(h(1e-2), r.subspace), kwargs...)
    prob = ODEProblem(shordinger(h, r.subspace; cache=cache), vec(r.state), (zero(t), t), h; save_everystep=false, save_start=false, alias_u0=true, kwargs...)
    result = solve(prob, algo)
    return r
end

end
