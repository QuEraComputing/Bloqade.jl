"""
    rydberg_density(reg, i::Int) -> Real

Calculates the rydberg density at site `i`.

```math
\\langle n_i \\rangle
```
"""
rydberg_density(reg, i::Int) = real(expect(put(nqubits(reg), i => Op.n), reg))

"""
    rydberg_density(reg) -> Vector

Return the rydberg density at each site.
"""
rydberg_density(reg) = [rydberg_density(reg, i) for i in 1:nqubits(reg)]

"""
    rydberg_corr([op=Op.n], reg) -> Matrix

Calculates the rydberg correlation matrix.

```math
\\langle \\text{op}_i \\text{op}_j \\rangle
```

here `op` can be `Op.n`, `X` or `Y`.

# Arguments

- `op`: the correlation function, default is `Op.n`.
- `reg`: required, the register object.
"""
rydberg_corr(reg) = rydberg_corr(Op.n, reg)

function rydberg_corr(op, reg)
    return [expect(chain(nqubits(reg), put(i => op), put(j => op)), reg) for i in 1:nqubits(reg), j in 1:nqubits(reg)]
end


"""
    get_average_rydberg_densities(atoms, reg; [C=2π * 862690 * MHz*µm^6], Ω[, ϕ, Δ], [dt=1e-3 * μs])

Return average Rydberg densities throughout an evolution

# Arguments

- `atoms`: a collection of atom positions.
- `reg`: required, the register object.

# Keyword Arguments

- `C`: optional, default unit is `MHz*µm^6`, interation parameter,
    see also [`RydInteract`](@ref).
- `Ω`: optional, default unit is `MHz`, Rabi frequencies, divided by 2, see also [`SumOfX`](@ref).
- `ϕ`: optional, does not have unit, the phase, see [`SumOfXPhase`](@ref).
- `Δ`: optional, default unit is `MHz`, detuning parameter, see [`SumOfN`](@ref).
- `dt`: optional, default unit is `μs`, time step for the evolution
"""
function get_average_rydberg_densities(atoms, reg; C = 2π * 862690, Ω = nothing, ϕ = nothing, Δ = nothing, dt=1e-3)

    # Get the duration for the evolution
    if isnothing(Ω) && isnothing(ϕ) && isnothing(Δ)
        error("At least one of Ω, ϕ or Δ needs to be specified for determining the duration of the evolution.")
    end

    allwaveforms = [Ω, ϕ, Δ]
    allwaveforms = allwaveforms[allwaveforms.!=nothing]
    allwaveforms = reduce(vcat, allwaveforms)
    duration = allwaveforms[1].duration
    for i = 2 : length(allwaveforms)
        if allwaveforms[i].duration != duration
            error("The durations of waveforms are not consisteng.")
        end
    end

    # Start the evolution 
    h = rydberg_h(atoms, C=C, Δ=Δ, Ω=Ω, ϕ=ϕ)
    prob = SchrodingerProblem(reg, duration, h, progress=true);
    integrator = init(prob, Vern8());
    densities = []
    for _ in TimeChoiceIterator(integrator, 0.0:dt:duration)
        normalize!(reg) 
        push!(densities, rydberg_density(reg)) 
    end

    return densities    
end