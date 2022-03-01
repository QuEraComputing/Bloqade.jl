"""
    simple_rydberg(n::Int, ϕ::Number)

Create a simple rydberg hamiltonian that has only [`XTerm`](@ref).
"""
simple_rydberg(n::Int, ϕ::Number) = XTerm(n, one(ϕ), ϕ)


function rydberg_h(atoms, C, Ω, ϕ, Δ)
    return RydInteract(atoms, C) + XTerm(length(atoms), Ω, ϕ) - NTerm(length(atoms), Δ)
end

function rydberg_h(atoms, Ω, ϕ, Δ)
    return RydInteract(atoms) + XTerm(length(atoms), Ω, ϕ) - NTerm(length(atoms), Δ)
end

"""
    rydberg_h(atoms; [C=2π * 862690 * MHz*µm^6], Ω[, ϕ, Δ])

Create a rydberg hamiltonian

```math
∑ \\frac{C}{|r_i - r_j|^6} n_i n_j + \\frac{Ω}{2} σ_x - Δ σ_n
```

shorthand for

```julia
RydInteract(C, atoms) + XTerm(length(atoms), Ω, ϕ) + ZTerm(length(atoms), Δ)
```

# Arguments

- `atoms`: a collection of atom positions

# Keyword Arguments

- `C`: optional, default unit is `MHz*µm^6`, interation parameter,
    see also [`RydInteract`](@ref).
- `Ω`: required, default unit is `MHz`, Rabi frequencies, see [`XTerm`](@ref).
- `Δ`: optional, default unit is `MHz`, detuning parameter, see [`NTerm`](@ref).
- `ϕ`: optional, does not have unit, the phase, see [`XTerm`](@ref).

!!! tips

    The parameters of Hamiltonian have their own default units to match hardware,
    one can use [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)
    to specify their units explicitly. If the units are specified explicitly,
    they will be converted to default units automatically.

# Example

```julia-repl
julia> using EaRyd

julia> atoms = generate_sites(SquareLattice(), 3, 3);

julia> rydberg_h(atoms; Δ=1.2, Ω=1.1)
Hamiltonian
  Term 1
   ∑(n=1:9) 1.1/2 σ^x

  Term 2
   -∑(n=1:9) 1.2 n

  Term 3
   ∑(n=1:9) 686.0/|r_i - r_j|^6 n_i n_j
```

```julia-repl
julia> rydberg_h(atoms; Δ=1.2, Ω=1.1, ϕ=2.1)
Hamiltonian
  Term 1
   ∑(n=1:9) 1.1/2 (e^{2.1i}|0)⟨1| + e^{-2.1i}|1⟩⟨0|)

  Term 2
   -∑(n=1:9) 1.2 n

  Term 3
   ∑(n=1:9) 686.0/|r_i - r_j|^6 n_i n_j
```
"""
function rydberg_h(atoms; Ω, Δ=nothing, C=nothing, ϕ=nothing) 
    if eltype(atoms) <: Tuple
        atoms = map(atoms) do each
            RydAtom(each)
        end
        atoms = SVector{length(atoms)}(atoms)
    end

    if ϕ === nothing
        term = XTerm(length(atoms), Ω)
    else
        term = XTerm(length(atoms), Ω, ϕ)
    end

    if Δ !== nothing
        term -= NTerm(length(atoms), Δ)
    end

    if C === nothing
        term += RydInteract(atoms)
    else
        term += RydInteract(atoms, C)
    end
    return term
end

function is_time_dependent(h::XTerm)
    return !(h.Ωs isa ConstParamType || isnothing(h.Ωs)) ||
        !(h.ϕs isa ConstParamType || isnothing(h.ϕs))
end

function is_time_dependent(h::Union{ZTerm, NTerm})
    return !(h.Δs isa ConstParamType)
end

# NOTE: we currently assume atom positions are constant
# might change in the future

is_time_dependent(h::RydInteract) = false
is_time_dependent(t::Negative) = is_time_dependent(t.term)

"""
    is_time_dependent(h::Hamiltonian)

Check if a hamiltonian is time-dependent.
"""
function is_time_dependent(h::Hamiltonian)
    any(is_time_dependent, h.terms)
end

# Time dependent Term
attime(x::Nothing, t::Real) = nothing
attime(x::Number, t::Real) = x
attime(x, t::Real) = x(t)
attime(x::AbstractArray, t::Real) = attime.(x, t)
attime(x::Tuple, t::Real) = attime.(x, t)

function (tm::XTerm)(t::Real)
    return XTerm(tm.nsites, attime(tm.Ωs, t), attime(tm.ϕs, t))
end

function (tm::ZTerm)(t::Real)
    return ZTerm(tm.nsites, attime(tm.Δs, t))
end

function (tm::NTerm)(t::Real)
    return NTerm(tm.nsites, attime(tm.Δs, t))
end

(tm::Negative)(t::Real) = Negative(tm.term(t))

function (tm::Hamiltonian)(t::Real)
    return Hamiltonian(map(x->x(t), tm.terms))
end

# fallback to constants
(tm::AbstractTerm)(t::Real) = tm
