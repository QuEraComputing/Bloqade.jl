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
function rydberg_h(atom_positions; C = 2π * 862690, Ω = nothing, ϕ = nothing, Δ = nothing)
    return rydberg_h(atom_positions, C, Ω, ϕ, Δ)
end

function rydberg_h(atom_positions, C, Ω, ϕ, Δ)
    positions = map(atom_positions) do pos
        return (pos...,)
    end

    nsites = length(positions)
    rydberg_term = RydInteract(positions, C)

    Ω = div_by_two(Ω)
 
    if !isnothing(Ω) && !isnothing(ϕ)
        rabi_term = SumOfXPhase(nsites, Ω, ϕ)
    elseif !isnothing(Ω) && isnothing(ϕ)
        rabi_term = SumOfX(nsites, Ω)
    elseif isnothing(Ω) && !isnothing(ϕ)
        @warn "Rydberg Hamiltonian contains non-zero rabi phase ϕ with no rabi amplitude Ω."
        rabi_term = SumOfXPhase(nsites, 0, ϕ)
    else
        rabi_term = nothing
    end

    if !isnothing(Δ)
        detuning_term = SumOfN(nsites, Δ)
    else
        detuning_term = nothing
    end

    return RydbergHamiltonian(rydberg_term,rabi_term,detuning_term)
end

function rydberg_h_3(atom_positions; 
        C = 2π * 862690, Ω_hf = nothing, ϕ_hf = nothing, Δ_hf = nothing, 
        Ω_r = nothing, ϕ_r = nothing, Δ_r = nothing)
    return rydberg_h_3(atom_positions, C, Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r)
end

function rydberg_h_3(atom_positions, C, Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r)
    positions = map(atom_positions) do pos
        return (pos...,)
    end

    nsites = length(positions)
    rydberg_term = RydInteract(positions, C; nlevel = 3)

    Ω_hf = div_by_two(Ω_hf)
    if !isnothing(Ω_hf) && !isnothing(ϕ_hf)
        rabi_term_hf = SumOfXPhase(nsites, Ω_hf, ϕ_hf; nlevel = 3, name = :hyperfine)
    elseif !isnothing(Ω_hf) && isnothing(ϕ_hf)
        rabi_term_hf = SumOfX(nsites, Ω_hf; nlevel = 3, name = :hyperfine)
    elseif isnothing(Ω_hf) && !isnothing(ϕ_hf)
        @warn "Rydberg Hamiltonian contains non-zero rabi phase ϕ_hf with no rabi amplitude Ω_hf."
        rabi_term_hf = SumOfXPhase(nsites, 0, ϕ_hf; nlevel = 3, name = :hyperfine)
    else
        rabi_term_hf = nothing
    end
    
    Ω_r = div_by_two(Ω_r)
    if !isnothing(Ω_r) && !isnothing(ϕ_r)
        rabi_term_r =+ SumOfXPhase(nsites, Ω_r, ϕ_r; nlevel = 3, name = :rydberg)
    elseif !isnothing(Ω_r) && isnothing(ϕ_r)
        rabi_term_r =+ SumOfX(nsites, Ω_r; nlevel = 3, name = :rydberg)
    elseif isnothing(Ω_r) && !isnothing(ϕ_r)
        @warn "Rydberg Hamiltonian contains non-zero rabi phase ϕ_r with no rabi amplitude Ω_r."
        rabi_term_r =+ SumOfXPhase(nsites, 0, ϕ_r; nlevel = 3, name = :rydberg)
    else
        rabi_term_r = nothing
    end

    if isnothing(rabi_term_hf) 
        rabi_term = rabi_term_r
    elseif isnothing(rabi_term_r)
        rabi_term = rabi_term_hf
    else
        rabi_term = rabi_term_hf + rabi_term_r
    end

    if !isnothing(Δ_hf)
        detuning_term = SumOfN(nsites, Δ_hf; nlevel = 3, name = :hyperfine)
        if !isnothing(Δ_r)
            detuning_term += SumOfN(nsites, add_params(Δ_hf, Δ_r); nlevel = 3, name = :rydberg)
        else
            detuning_term += SumOfN(nsites, Δ_hf; nlevel = 3, name = :rydberg)
        end
    elseif !isnothing(Δ_r)
        detuning_term = SumOfN(nsites, Δ_r; nlevel = 3, name = :rydberg)
    else
        detuning_term = nothing
    end

    !isnothing(detuning_term) && (detuning_term = -detuning_term)

    rh3 = rydberg_term
    for t in (rabi_term, detuning_term)
        !isnothing(t) && (rh3 += t) 
    end

    return rh3
end

function add_params(a, b)
    if a isa Vector && b isa Vector
        return add_param.(a, b)
    elseif a isa Vector && !(b isa Vector)
        return add_param.(a, fill(b, length(a)))
    elseif !(a isa Vector) && b isa Vector
        return add_param.(fill(a, length(b)), b)
    else
        return add_param(a, b)
    end
end

function add_param(a, b)
    if is_const_param(a) && is_const_param(b)
        return a + b
    elseif is_time_function(a) && is_const_param(b)
        return x -> a(x) + b
    elseif is_const_param(a) && is_time_function(b)
        return x -> a + b(x)
    else
        return x -> a(x) + b(x)
    end
end

function div_by_two(Ω)
    isnothing(Ω) && return
    if is_const_param(Ω)
        return Ω ./ 2
    end

    return if Ω isa Vector
        map(Ω) do Ω_i
            return DivByTwo(Ω_i)
        end
    else
        DivByTwo(Ω)
    end
end

attime(t::Real) = h -> attime(h, t)

function attime(h::AbstractBlock, t::Real)
    blks = map(subblocks(h)) do blk
        return attime(blk, t)
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
function attime(h::SumOfX{D, name}, t::Real) where {D, name}
    is_const_param(h.Ω) && return h
    if h.Ω isa Vector
        SumOfX(h.nsites, map(x -> x(t), h.Ω); nlevel = D, name = name)
    else
        SumOfX(h.nsites, h.Ω(t); nlevel = D, name = name)
    end
end

function attime(h::SumOfXPhase{D, name}, t::Real) where {D, name}
    if is_const_param(h.Ω)
        Ω = h.Ω
    elseif h.Ω isa Vector
        Ω = map(x -> x(t), h.Ω)
    else
        Ω = h.Ω(t)
    end

    if is_const_param(h.ϕ)
        ϕ = h.ϕ
    elseif h.ϕ isa Vector
        ϕ = map(x -> x(t), h.ϕ)
    else
        ϕ = h.ϕ(t)
    end

    return SumOfXPhase(h.nsites, Ω, ϕ; nlevel = D, name = name)
end

function attime(h::Union{SumOfZ,SumOfN}, t::Real)
    if is_const_param(h.Δ)
        Δ = h.Δ
    elseif h.Δ isa Vector
        Δ = map(x -> x(t), h.Δ)
    else
        Δ = h.Δ(t)
    end
    return typeof(h)(h.nsites, Δ)
end


function attime(h::RydbergHamiltonian, t::Real)
    return attime(add_terms(h),t)
end