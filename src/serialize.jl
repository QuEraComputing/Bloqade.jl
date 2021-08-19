Configurations.is_option(::AbstractTerm) = true
Configurations.is_option(::Type{<:AbstractTerm}) = true

Configurations.alias(::Type{<:XTerm}) = "xterm"
Configurations.alias(::Type{<:ZTerm}) = "zterm"
Configurations.alias(::Type{<:NTerm}) = "nterm"
Configurations.alias(::Type{<:RydInteract}) = "rydberg"

function Configurations.field_default(::Type{<:AbstractTerm}, ::Symbol)
    return Configurations.no_default
end

function Configurations.to_dict(term::XTerm)
    d = OrderedDict{String, Any}(
        "nsites" => term.nsites,
    )

    isnothing(term.Ωs) || (d["omega"] = term.Ωs)
    isnothing(term.ϕs) || (d["phi"] = term.ϕs)
    return d
end

function Configurations.to_dict(term::ZTerm)
    d = OrderedDict{String, Any}(
        "nsites" => term.nsites,
    )
    isnothing(term.Δs) || (d["delta"] = term.Δs)
    return d
end

function Configurations.to_dict(term::NTerm)
    d = OrderedDict{String, Any}(
        "nsites" => term.nsites,
    )
    isnothing(term.Δs) || (d["delta"] = term.Δs)
    return d
end

function Configurations.to_dict(term::RydInteract)
    return OrderedDict{String, Any}(
        "atoms" => [collect(x) for x in term.atoms],
        "C" => term.C,
    )
end

function Configurations.to_dict(term::Hamiltonian)
    d = OrderedDict{String, Any}()
    for t in term.terms
        d[Configurations.alias(typeof(t))] = Configurations.to_dict(t)
    end
    return d
end

function validate_term_keys(d, available_keys)
    for k in keys(d)
        k == "#filename#" && continue
        k in available_keys || error("invalid key: $k")
    end
end

function Configurations.from_dict_validate(::Type{T}, d::AbstractDict{String}) where {T <: XTerm}
    validate_term_keys(d, ["Ωs", "ϕs", "omega", "phi", "nsites"])

    haskey(d, "Ωs") || haskey(d, "omega") || error("key Ωs/omega (Rabi frequency) is required")
    rabi = haskey(d, "Ωs") ? d["Ωs"] : d["omega"]
    
    if rabi isa Number
        haskey(d, "nsites") || error("key nsites is required for scalar Ωs/omega (Rabi frequency)")
    end

    # check duplicated keys
    haskey(d, "Ωs") && haskey(d, "omega") && error("key Ωs/omega (Rabi frequency) is duplicated")
    haskey(d, "ϕs") && haskey(d, "phi") && error("key ϕs/phi (phase) is duplicated")

    return Configurations.from_dict_inner(T, d)
end

function Configurations.from_dict_validate(::Type{T}, d::AbstractDict{String}) where {T <: ZTerm}
    validate_term_keys(d, ["Δs", "delta", "nsites"])
    haskey(d, "Δs") || haskey(d, "delta") || error("key Δs/delta (detuning) is required")
    detuning = haskey(d, "Δs") ? d["Δs"] : d["delta"]

    if detuning isa Number
        haskey(d, "nsites") || error("key nsites is required for scalar Δs/delta (detuning)")
    end

    # check duplicated keys
    haskey(d, "Δs") && haskey(d, "delta") && error("key Δs/delta (detuning) is duplicated")
    return Configurations.from_dict_inner(T, d)
end

function Configurations.from_dict_validate(::Type{T}, d::AbstractDict{String}) where {T <: NTerm}
    validate_term_keys(d, ["Δs", "delta", "nsites"])
    haskey(d, "Δs") || haskey(d, "delta") || error("key Δs/delta (detuning) is required")
    detuning = haskey(d, "Δs") ? d["Δs"] : d["delta"]

    if detuning isa Number
        haskey(d, "nsites") || error("key nsites is required for scalar Δs/delta (detuning)")
    end

    # check duplicated keys
    haskey(d, "Δs") && haskey(d, "delta") && error("key Δs/delta (detuning) is duplicated")
    return Configurations.from_dict_inner(T, d)
end

function Configurations.from_dict_validate(::Type{T}, d::AbstractDict{String}) where {T <: RydInteract}
    validate_term_keys(d, ["atoms", "C"])

    haskey(d, "atoms") || error("key atoms is required")
    haskey(d, "C") || error("key C is required")

    d["C"] isa Number || error("key C must be a Number")
    return Configurations.from_dict_inner(T, d)
end

function Configurations.from_dict_validate(::Type{T}, d::AbstractDict{String}) where {T <: Hamiltonian}
    validate_term_keys(d, ["xterm", "zterm", "nterm", "rydberg"])
    return Configurations.from_dict_inner(T, d)
end

function Configurations.from_dict_inner(::Type{T}, d::AbstractDict{String}) where {T <: XTerm}
    args = []
    haskey(d, "nsites") && push!(args, d["nsites"])
    haskey(d, "omega") && push!(args, d["omega"])
    haskey(d, "Ωs") && push!(args, d["Ωs"])
    haskey(d, "phi") && push!(args, d["phi"])
    haskey(d, "ϕs") && push!(args, d["ϕs"])
    return XTerm(args...)
end

function Configurations.from_dict_inner(::Type{T}, d::AbstractDict{String}) where {T <: ZTerm}
    args = []
    haskey(d, "nsites") && push!(args, d["nsites"])
    haskey(d, "Δs") && push!(args, d["Δs"])
    haskey(d, "delta") && push!(args, d["delta"])
    return ZTerm(args...)
end

function Configurations.from_dict_inner(::Type{T}, d::AbstractDict{String}) where {T <: NTerm}
    args = []
    haskey(d, "nsites") && push!(args, d["nsites"])
    haskey(d, "Δs") && push!(args, d["Δs"])
    haskey(d, "delta") && push!(args, d["delta"])
    return NTerm(args...)
end

function Configurations.from_dict_inner(::Type{T}, d::AbstractDict{String}) where {T <: RydInteract}
    if d["atoms"] isa String
        path = haskey(d, "#filename#") ? joinpath(d["#filename#"], d["atoms"]) : d["atoms"]
        atoms = read_atoms(path)
    elseif d["atoms"] isa Vector
        atoms = [RydAtom(x) for x in d["atoms"]]
    else
        error("invalid type for field atoms")
    end
    return RydInteract(atoms, d["C"])
end

function Configurations.from_dict_inner(::Type{T}, d::AbstractDict{String}) where {T <: Hamiltonian}
    terms = []
    haskey(d, "xterm") && push!(terms, from_dict(XTerm, d["xterm"]))
    haskey(d, "zterm") && push!(terms, from_dict(ZTerm, d["zterm"]))
    haskey(d, "nterm") && push!(terms, from_dict(NTerm, d["nterm"]))
    haskey(d, "rydberg") && push!(terms, from_dict(RydInteract, d["rydberg"]))
    return Hamiltonian((terms..., ))
end

"""
    write_atoms(io::IO, atoms::AbstractVector{<:RydAtom})

Write a list of atom positions to stream `io`.
"""
function write_atoms(io::IO, atoms::AbstractVector{<:RydAtom})
    for atom in atoms
        println(io, join(map(x->string(x), atom), "  "))
    end
    return
end

"""
    write_atoms(filename::String, atoms::AbstractVector{<:RydAtom})

Write a list of atom positions to the file given by `filename`.

# Example

This saves the atom position generated from [`square_lattice`](@ref)
to a file `demo.atoms`.

```julia
atoms = square_lattice(5, 0.8)
write_atoms("demo.atoms", atoms)
```
"""
function write_atoms(filename::String, atoms::AbstractVector{<:RydAtom})
    return open(filename, "w+") do io
        write_atoms(io, atoms)
    end
end

read_atoms(io::IO) = read_atoms(io, Int)

"""
    read_atoms(io::IO[, T=Int])

Read atom positions from stream `io`.
"""
function read_atoms(io::IO, ::Type{T}) where T
    atoms = [RydAtom(row) for row in eachrow(readdlm(io, T))]
    return SVector{length(atoms)}(atoms)
end

"""
    read_atoms(filename[, T=Int])

Read atom positions from file `filename`.
"""
function read_atoms(filename::String, T...)
    return open(filename) do io
        read_atoms(io, T...)
    end
end
