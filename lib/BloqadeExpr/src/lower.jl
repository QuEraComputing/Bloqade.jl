function emit_dynamic_terms(ex::Add)
    list = []
    for t in subblocks(ex)
        append!(list, emit_dynamic_terms(t))
    end
    return list
end

function emit_dynamic_terms(ex::Scale)
    return map(emit_dynamic_terms(content(ex))) do (f, t)
        return f => factor(ex) * t
    end
end

emit_dynamic_terms(ex::RydInteract) = Any[one=>ex,]

function emit_dynamic_terms(ex::SumOfX)
    return if ex.Ω isa Number || ex.Ω isa Vector{<:Number}
        Any[one=>ex,]
    elseif ex.Ω isa Vector && is_time_function(ex.Ω)
        Any[Ω => put(ex.nsites, i => X) for (i, Ω) in enumerate(ex.Ω)]
    elseif is_time_function(ex.Ω)
        Any[ex.Ω=>sum(put(ex.nsites, i => X) for i in 1:ex.nsites)]
    else
        error("unexpected type for Ω: $(typeof(ex.Ω))")
    end
end

function emit_dynamic_terms(ex::Union{SumOfN,SumOfZ})
    op = ex isa SumOfN ? ConstGate.P1 : Z

    return if ex.Δ isa Number || ex.Δ isa Vector{<:Number}
        Any[one=>ex]
    elseif ex.Δ isa Vector && is_time_function(ex.Δ)
        Any[Δ => put(ex.nsites, i => op) for (i, Δ) in enumerate(ex.Δ)]
    elseif is_time_function(ex.Δ)
        Any[ex.Δ=>sum(put(ex.nsites, i => op) for i in 1:ex.nsites)]
    else
        error("unexpected type for Δ: $(typeof(ex.Δ))")
    end
end

function emit_dynamic_terms(ex::SumOfXPhase)
    @switch (ex.Ω, ex.ϕ) begin
        @case (::Number, ::Number) || (::Vector{<:Number}, ::Vector{<:Number})
        return Any[one=>ex]

        # Ω time-dependent
        @case (Ω::Vector, ϕ::Number)
        return map(enumerate(Ω)) do (i, Ω_i)
            return Ω_i => put(ex.nsites, i => XPhase(ϕ))
        end
        @case (Ω, ϕ::Number)
        return Any[Ω=>SumOfXPhase(ex.nsites, one(ϕ), ϕ)]

        @case (Ω::Vector, ϕ::Vector{<:Number})
        return map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
            return Ω_i => put(ex.nsites, i => XPhase(ϕ_i))
        end
        @case (Ω, ϕ::Vector{<:Number})
        return Any[Ω=>SumOfXPhase(ex.nsites, one(eltype(ϕ)), ϕ)]

        # ϕ time-dependent
        @case (Ω::Number, ϕ::Vector)
        lhs = map(enumerate(ϕ)) do (i, ϕ_i)
            return (t -> exp(ϕ_i(t) * im)) => put(ex.nsites, i => (Ω * ConstGate.Pu))
        end

        rhs = map(enumerate(ϕ)) do (i, ϕ_i)
            return (t -> exp(-ϕ_i(t) * im)) => put(ex.nsites, i => (Ω * ConstGate.Pd))
        end
        return vcat(lhs, rhs)

        @case (Ω::Number, ϕ)
        return [
            (t -> exp(ϕ(t) * im)) => sum(put(ex.nsites, i => (Ω * ConstGate.Pu)) for i in 1:ex.nsites),
            (t -> exp(-ϕ(t) * im)) => sum(put(ex.nsites, i => (Ω * ConstGate.Pd)) for i in 1:ex.nsites),
        ]

        @case (Ω::Vector{<:Number}, ϕ::Vector)
        lhs = map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
            return (t -> exp(ϕ_i(t) * im)) => put(ex.nsites, i => (Ω_i * ConstGate.Pu))
        end

        rhs = map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
            return (t -> exp(-ϕ_i(t) * im)) => put(ex.nsites, i => (Ω_i * ConstGate.Pd))
        end
        return vcat(lhs, rhs)
        @case (Ω::Vector{<:Number}, ϕ)
        return [
            (t -> exp(ϕ(t) * im)) => sum(put(ex.nsites, i => (Ω_i * ConstGate.Pu)) for (i, Ω_i) in enumerate(Ω)),
            (t -> exp(-ϕ(t) * im)) => sum(put(ex.nsites, i => (Ω_i * ConstGate.Pd)) for (i, Ω_i) in enumerate(Ω)),
        ]

        # both time-dependent
        @case (Ω::Vector, ϕ::Vector)
        lhs = map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
            return (t -> Ω_i(t) * exp(ϕ_i(t) * im)) => put(ex.nsites, i => ConstGate.Pu)
        end

        rhs = map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
            return (t -> Ω_i(t) * exp(-ϕ_i(t) * im)) => put(ex.nsites, i => ConstGate.Pd)
        end
        return vcat(lhs, rhs)
        @case (Ω::Vector, ϕ)
        lhs = map(enumerate(Ω)) do (i, Ω_i)
            return (t -> Ω_i(t) * exp(ϕ(t) * im)) => put(ex.nsites, i => ConstGate.Pu)
        end

        rhs = map(enumerate(Ω)) do (i, Ω_i)
            return (t -> Ω_i(t) * exp(-ϕ(t) * im)) => put(ex.nsites, i => ConstGate.Pd)
        end
        return vcat(lhs, rhs)
        @case (Ω, ϕ::Vector)
        lhs = map(enumerate(ϕ)) do (i, ϕ_i)
            return (t -> Ω(t) * exp(ϕ_i(t) * im)) => put(ex.nsites, i => ConstGate.Pu)
        end

        rhs = map(enumerate(ϕ)) do (i, ϕ_i)
            return (t -> Ω(t) * exp(-ϕ_i(t) * im)) => put(ex.nsites, i => ConstGate.Pd)
        end
        return vcat(lhs, rhs)
        @case (Ω, ϕ)
        return [
            (t -> Ω(t) * exp(ϕ(t) * im)) => sum(put(ex.nsites, i => ConstGate.Pu) for i in 1:ex.nsites),
            (t -> Ω(t) * exp(-ϕ(t) * im)) => sum(put(ex.nsites, i => ConstGate.Pd) for i in 1:ex.nsites),
        ]
    end
end

"""
    Hamiltonian(::Type{Tv}, expr[, space=fullspace])

Create a `Hamiltonian` from hamiltonian expr that has matrix element of type `Tv`.
"""
function Hamiltonian(::Type{Tv}, ex::AbstractBlock, space::AbstractSpace = fullspace) where {Tv}
    fs, ts = [], []
    const_term = nothing
    for (f, op) in emit_dynamic_terms(ex)
        if f === Base.one
            const_term = isnothing(const_term) ? op : const_term + op
        else
            push!(fs, f)
            push!(ts, mat(Tv, op, space))
        end
    end
    push!(fs, Base.one)
    isnothing(const_term) || push!(ts, mat(Tv, const_term, space))
    return Hamiltonian((fs...,), (ts...,))
end

function YaoBlocks.Optimise.to_basictypes(h::SumOfX)
    is_const_param(h.Ω) || throw(ArgumentError("expect constant hamiltonian"))
    return if h.Ω isa Vector
        sum(h.Ω[i] * put(h.nsites, i => X) for i in 1:h.nsites)
    else
        h.Ω * sum(put(h.nsites, i => X) for i in 1:h.nsites)
    end
end

function YaoBlocks.Optimise.to_basictypes(h::Union{SumOfN,SumOfZ})
    is_const_param(h.Δ) || throw(ArgumentError("expect constant hamiltonian"))
    op = h isa SumOfN ? ConstGate.P1 : Z
    return if h.Δ isa Vector
        sum(h.Δ[i] * put(h.nsites, i => op) for i in 1:h.nsites)
    else
        h.Δ * sum(put(h.nsites, i => op) for i in 1:h.nsites)
    end
end

function YaoBlocks.Optimise.to_basictypes(h::SumOfXPhase)
    is_const_param(h.Ω) || throw(ArgumentError("expect constant hamiltonian"))
    is_const_param(h.ϕ) || throw(ArgumentError("expect constant hamiltonian"))
    return @switch (h.Ω, h.ϕ) begin
        @case (::Vector, ::Vector)
        sum(h.Ω[i] * put(h.nsites, i => XPhase(h.ϕ[i])) for i in 1:h.nsites)
        @case (Ω::Vector, ϕ)
        sum(Ω_i * put(h.nsites, i => XPhase(ϕ)) for (i, Ω_i) in enumerate(Ω))
        @case (Ω, ϕ::Vector)
        sum(Ω * put(h.nsites, i => XPhase(ϕ_i)) for (i, ϕ_i) in enumerate(ϕ))
        @case (Ω, ϕ)
        sum(Ω * put(h.nsites, i => XPhase(ϕ)) for i in 1:h.nsites)
    end
end

function YaoBlocks.Optimise.to_basictypes(ex::RydInteract)
    nsites = length(ex.atoms)

    term = nothing
    for i in 1:nsites, j in 1:i-1
        x, y = ex.atoms[i], ex.atoms[j]
        h = ex.C / distance(x, y)^6 * kron(nsites, i => ConstGate.P1, j => ConstGate.P1)

        if isnothing(term)
            term = h
        else
            term += h
        end
    end
    term === nothing && return Add(nsites)
    return term
end

function emit_lowered(h)
    return YaoBlocks.Optimise.simplify(h; rules = [YaoBlocks.Optimise.to_basictypes])
end
