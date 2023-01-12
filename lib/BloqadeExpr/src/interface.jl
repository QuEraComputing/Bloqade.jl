module Op

using YaoBlocks
const n = YaoBlocks.ConstGate.P1

end

process_atom_positions(atom_positions::BoundedLattice) = atom_positions

function process_atom_positions(atom_positions) # catches other cases
    return map(atom_positions) do pos
        (pos...,)
    end
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
∑ \\frac{C}{|x_i - x_j|^6} n_i n_j + \\frac{Ω}{2} σ_x - Δ σ_n
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
∑ 5.42e6/|x_i-x_j|^6 n_i n_j
```

```julia-repl
julia> rydberg_h(atoms; Ω=0.1)
nqubits: 4
+
├─ ∑ 5.42e6/|x_i-x_j|^6 n_i n_j
└─ 0.05 ⋅ ∑ σ^x_i
```
"""
function rydberg_h(atom_positions; C = 2π * 862690, Ω = nothing, ϕ = nothing, Δ = nothing)
    return rydberg_h(atom_positions, C, Ω, ϕ, Δ)
end

function rydberg_h(atom_positions, C, Ω, ϕ, Δ)
    positions = process_atom_positions(atom_positions)

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

"""
    rydberg_h_3(atoms; [C=2π * 862690 * MHz*µm^6, 
        Ω_hf = nothing, ϕ_hf = nothing, Δ_hf = nothing, 
        Ω_r = nothing, ϕ_r = nothing, Δ_r = nothing])

Create a 3-level Rydberg Hamiltonian

```math
\\sum_{i<j} \\frac{C}{|x_i - x_j|^6} n^r_i n^r_j + 
\\sum_{i} \\left[\\frac{Ω^{\\mathrm{hf}}}{2} (e^{iϕ^\\mathrm{hf}}|0⟩⟨1| + e^{-iϕ^\\mathrm{hf}}|1⟩⟨0|) - Δ^{\\mathrm{hf}} n^{1}_i + 
\\frac{Ω^{\\mathrm{r}}}{2} (e^{iϕ^\\mathrm{r}}|1⟩⟨r| + e^{-iϕ^\\mathrm{r}}|r⟩⟨1|) - (Δ^{\\mathrm{hf}} + Δ^{\\mathrm{r}}) n^{\\mathrm{r}}_i \\right]
```

shorthand for

```julia
RydInteract(C, atoms; nlevel = 3) + 
    SumOfXPhase_01(length(atoms), Ω_hf/2, ϕ_hf) - SumOfN(length(atoms), Δ_hf) +
    SumOfXPhase_1r(length(atoms), Ω_r/2, ϕ_r) - SumOfN(length(atoms), Δ_r + Δ_hf)
```
"""
function rydberg_h_3(atom_positions; 
        C = 2π * 862690, Ω_hf = nothing, ϕ_hf = nothing, Δ_hf = nothing, 
        Ω_r = nothing, ϕ_r = nothing, Δ_r = nothing)
    return rydberg_h_3(atom_positions, C, Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r)
end

function rydberg_h_3(atom_positions, C, Ω_hf, ϕ_hf, Δ_hf, Ω_r, ϕ_r, Δ_r)
    positions = process_atom_positions(atom_positions)

    nsites = length(positions)
    rydberg_term = RydInteract(positions, C; nlevel = 3)

    Ω_hf = div_by_two(Ω_hf)
    if !isnothing(Ω_hf) && !isnothing(ϕ_hf)
        rabi_term_hf = SumOfXPhase_01(nsites, Ω_hf, ϕ_hf)
    elseif !isnothing(Ω_hf) && isnothing(ϕ_hf)
        rabi_term_hf = SumOfX_01(nsites, Ω_hf)
    elseif isnothing(Ω_hf) && !isnothing(ϕ_hf)
        @warn "Rydberg Hamiltonian contains non-zero rabi phase ϕ_hf with no rabi amplitude Ω_hf."
        rabi_term_hf = SumOfXPhase_01(nsites, 0, ϕ_hf)
    else
        rabi_term_hf = nothing
    end
    
    Ω_r = div_by_two(Ω_r)
    if !isnothing(Ω_r) && !isnothing(ϕ_r)
        rabi_term_r = SumOfXPhase_1r(nsites, Ω_r, ϕ_r)
    elseif !isnothing(Ω_r) && isnothing(ϕ_r)
        rabi_term_r = SumOfX_1r(nsites, Ω_r)
    elseif isnothing(Ω_r) && !isnothing(ϕ_r)
        @warn "Rydberg Hamiltonian contains non-zero rabi phase ϕ_r with no rabi amplitude Ω_r."
        rabi_term_r = SumOfXPhase_1r(nsites, 0, ϕ_r)
    else
        rabi_term_r = nothing
    end

    detuning_term_hf = !isnothing(Δ_hf) ? SumOfN_1(nsites, Δ_hf) : nothing
    detuning_term_r = if !isnothing(Δ_r)
        if !isnothing(Δ_hf) 
            SumOfN_r(nsites, add_params(Δ_hf, Δ_r))
        else
            SumOfN_r(nsites, Δ_r)
        end
    elseif !isnothing(Δ_hf)
        SumOfN_r(nsites, Δ_hf)
    else
        nothing
    end

    rh3 = RydbergHamiltonian3(rydberg_term, rabi_term_hf, detuning_term_hf, rabi_term_r, detuning_term_r)
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



function mult_by_two(Ω)
    isnothing(Ω) && return
    if is_const_param(Ω)
        return Ω .* 2
    end

    return if Ω isa Vector
        map(Ω) do Ω_i
            return Ω_i.f
        end
    else
        Ω.f
    end

end

"""
    get_rydberg_params(h)

Returns a named tuple containing the fields `atoms`, `ϕ`, `Ω`, and `Δ` for the RydbergHamiltonian, `h`.

See also [`rydberg_h`](@ref)

# Example

```julia
    (atoms,ϕ,Ω,Δ) = get_rydberg_params(h)
````

"""
function get_rydberg_params(h::RydbergHamiltonian)
    # extracts parameters from RydbergHamiltonian
    ϕ = nothing
    Ω = nothing
    Δ = nothing

    if h.rabi_term isa SumOfX
        Ω = mult_by_two(h.rabi_term.Ω)
    elseif h.rabi_term isa SumOfXPhase
        Ω = mult_by_two(h.rabi_term.Ω)
        ϕ = h.rabi_term.ϕ
    end

    if h.detuning_term isa SumOfN
        Δ = h.detuning_term.Δ
    end

    return (atoms=h.rydberg_term.atoms,ϕ=ϕ,Ω=Ω,Δ=Δ)
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
├─ ∑ 5.42e6/|x_i-x_j|^6 n_i n_j
└─ Ω(t) ⋅ ∑ σ^x_i

julia> h |> attime(0.1)
nqudits: 4
+
├─ ∑ 5.42e6/|x_i-x_j|^6 n_i n_j
└─ 0.0499 ⋅ ∑ σ^x_i
```
"""
attime(h::PrimitiveBlock, ::Real) = h
function attime(h::SumOfXTypes, t::Real)
    is_const_param(h.Ω) && return h
    T = typeof(h)
    if h.Ω isa Vector
        T(h.nsites, map(x -> x(t), h.Ω))
    else
        T(h.nsites, h.Ω(t))
    end
end

function attime(h::SumOfXPhaseTypes, t::Real)
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

    T = typeof(h)
    return T(h.nsites, Ω, ϕ)
end

function attime(h::SumOfZAndNTypes, t::Real)
    if is_const_param(h.Δ)
        Δ = h.Δ
    elseif h.Δ isa Vector
        Δ = map(x -> x(t), h.Δ)
    else
        Δ = h.Δ(t)
    end
    return typeof(h)(h.nsites, Δ)
end

attime(h::RydbergHamiltonian, t::Real) = attime(add_terms(h),t)
attime(h::RydbergHamiltonian3, t::Real) = attime(add_terms(h),t)
