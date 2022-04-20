module Op

using YaoBlocks
const n = YaoBlocks.ConstGate.P1

end

"""
    emulate!(prob)

Run emulation of a given problem.
"""
function emulate! end

"""
    rydberg_h(atoms; [C=2π * 862690 * MHz*µm^6], Ω[, ϕ, Δ])

Create a rydberg hamiltonian

```math
∑ \\frac{C}{|r_i - r_j|^6} n_i n_j + \\frac{Ω}{2} σ_x - Δ σ_n
```

shorthand for

```julia
RydInteract(C, atoms) + SumOfXPhase(length(atoms), Ω, ϕ) - SumOfN(length(atoms), Δ)
```

# Arguments

- `atoms`: a collection of atom positions.

# Keyword Arguments

- `C`: optional, default unit is `MHz*µm^6`, interation parameter,
    see also [`RydInteract`](@ref).
- `Ω`: optional, default unit is `MHz`, Rabi frequencies, divided by 2, see also [`SumOfX`](@ref).
- `Δ`: optional, default unit is `MHz`, detuning parameter, see [`SumOfN`](@ref).
- `ϕ`: optional, does not have unit, the phase, see [`SumOfXPhase`](@ref).

!!! tips
    The rabi frequencies are divided by two in the Rydberg hamiltonian unlike
    directly constructing via [`SumOfX`](@ref) or [`SumOfXPhase`](@ref).

!!! tips

    The parameters of Hamiltonian have their own default units to match hardware,
    one can use [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)
    to specify their units explicitly. If the units are specified explicitly,
    they will be converted to default units automatically.

# Example

```julia-repl
julia> using Bloqade

julia> atoms = [(1, ), (2, ), (3, ), (4, )]
4-element Vector{Tuple{Int64}}:
 (1,)
 (2,)
 (3,)
 (4,)

julia> rydberg_h(atoms)
∑ 5.42e6/|r_i-r_j|^6 n_i n_j
```

```julia-repl
julia> rydberg_h(atoms; Ω=0.1)
nqubits: 4
+
├─ ∑ 5.42e6/|r_i-r_j|^6 n_i n_j
└─ 0.05 ⋅ ∑ σ^x_i
```
"""
function rydberg_h(atom_positions; C=2π * 862690, Ω=nothing, ϕ=nothing, Δ=nothing)
    positions = map(atom_positions) do pos
        (pos..., )
    end

    nsites = length(positions)
    term = RydInteract(positions, C)

    Ω = div_by_two(Ω)

    if !isnothing(Ω) && !isnothing(ϕ)
        term += SumOfXPhase(nsites, Ω, ϕ)
    elseif !isnothing(Ω) && isnothing(ϕ)
        term += SumOfX(nsites, Ω)
    end

    if !isnothing(Δ)
        term -= SumOfN(nsites, Δ)
    end

    return YaoBlocks.Optimise.simplify(term)
end

function div_by_two(Ω)
    isnothing(Ω) && return
    if is_const_param(Ω)
        return Ω ./ 2
    end

    return if Ω isa Vector
        map(Ω) do Ω_i
            t->Ω_i(t) / 2
        end
    else
        t->Ω(t)/2
    end
end

attime(t::Real) = h->attime(h, t)

function attime(h::AbstractBlock, t::Real)
    blks = map(subblocks(h)) do blk
        attime(blk, t)
    end
    return chsubblocks(h, blks)
end

# Yao Blocks cannot take time-dependent function

"""
    attime(h, t)
    attime(t)

Return the hamiltonian at time `t`.

# Example

```julia
julia> h = rydberg_h(atoms; Ω=sin)
nqudits: 4
+
├─ ∑ 5.42e6/|r_i-r_j|^6 n_i n_j
└─ Ω(t) ⋅ ∑ σ^x_i

julia> h |> attime(0.1)
nqudits: 4
+
├─ ∑ 5.42e6/|r_i-r_j|^6 n_i n_j
└─ 0.0499 ⋅ ∑ σ^x_i
```
"""
attime(h::PrimitiveBlock, ::Real) = h
function attime(h::SumOfX, t::Real)
    is_const_param(h.Ω) && return h
    if h.Ω isa Vector
        SumOfX(h.nsites, map(x->x(t), h.Ω))
    else
        SumOfX(h.nsites, h.Ω(t))
    end
end

function attime(h::SumOfXPhase, t::Real)
    if is_const_param(h.Ω)
        Ω = h.Ω
    elseif h.Ω isa Vector
        Ω = map(x->x(t), h.Ω)
    else
        Ω = h.Ω(t)
    end

    if is_const_param(h.ϕ)
        ϕ = h.ϕ
    elseif h.ϕ isa Vector
        ϕ = map(x->x(t), h.ϕ)
    else
        ϕ = h.ϕ(t)
    end

    return SumOfXPhase(h.nsites, Ω, ϕ)
end

function attime(h::Union{SumOfZ, SumOfN}, t::Real)
    if is_const_param(h.Δ)
        Δ = h.Δ
    elseif h.Δ isa Vector
        Δ = map(x->x(t), h.Δ)
    else
        Δ = h.Δ(t)
    end
    return typeof(h)(h.nsites, Δ)
end
