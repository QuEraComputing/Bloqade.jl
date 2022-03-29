function emit_dynamic_terms(ex::Add)
    list = []
    for t in subblocks(ex)
        append!(list, emit_dynamic_terms(t))
    end
    return list
end

function emit_dynamic_terms(ex::Scale)
    return map(emit_dynamic_terms(content(ex))) do (f, t)
        f => factor(ex) * t
    end
end

emit_dynamic_terms(ex::RydInteract) = Any[one=>ex, ]

function emit_dynamic_terms(ex::SumOfX)
    return if ex.Ω isa Number || ex.Ω isa Vector{<:Number}
        Any[one=>ex, ]
    elseif ex.Ω isa Vector && is_time_function(ex.Ω)
        Any[Ω=>put(ex.nsites, i=>X) for (i, Ω) in enumerate(ex.Ω)]
    elseif is_time_function(ex.Ω)
        Any[ex.Ω=>sum(put(ex.nsites, i=>X) for i in 1:ex.nsites)]
    else
        error("unexpected type for Ω: $(typeof(ex.Ω))")
    end
end

function emit_dynamic_terms(ex::Union{SumOfN, SumOfZ})
    op = ex isa SumOfN ? ConstGate.P1 : Z

    return if ex.Δ isa Number || ex.Δ isa Vector{<:Number}
        Any[one=>ex]
    elseif ex.Δ isa Vector && is_time_function(ex.Δ)
        Any[Δ=>put(ex.nsites, i=>op) for (i, Δ) in enumerate(ex.Δ)]
    elseif is_time_function(ex.Δ)
        Any[ex.Δ=>sum(put(ex.nsites, i=>op) for i in 1:ex.nsites)]
    else
        error("unexpected type for Δ: $(typeof(ex.Δ))")
    end
end

function emit_dynamic_terms(ex::SumOfXPhase)
    @switch (ex.Ω, ex.ϕ) begin
        @case (::Number, ::Number) || (::Vector{<:Number}, ::Vector{<:Number})
            return Any[one=>ex]

        # Ω time-dependent
        @case (Ω::Vector, ϕ::Number) && if is_time_function(Ω) end
            return map(enumerate(Ω)) do (i, Ω_i)
                Ω_i => put(ex.nsites, i=>XPhase(ϕ))
            end
        @case (Ω, ϕ::Number) && if is_time_function(Ω) end
            return Any[Ω=>SumOfXPhase(ex.nsites, one(ϕ), ϕ)]

        @case (Ω::Vector, ϕ::Vector{<:Number}) && if is_time_function(Ω) end
            return map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
                Ω_i => put(ex.nsites, i=>XPhase(ϕ_i))
            end
        @case (Ω, ϕ::Vector{<:Number}) && if is_time_function(Ω) end
            return Any[Ω=>SumOfXPhase(ex.nsites, one(eltype(ϕ)), ϕ)]

        # ϕ time-dependent
        @case (Ω::Number, ϕ::Vector) && if is_time_function(ϕ) end
            lhs = map(enumerate(ϕ)) do (i, ϕ_i)
                (t->exp(ϕ_i(t) * im)) => put(ex.nsites, i=>matblock([0 Ω;0 0]))
            end

            rhs = map(enumerate(ϕ)) do (i, ϕ_i)
                (t->exp(-ϕ_i(t) * im)) => put(ex.nsites, i=>matblock([0 0;Ω 0]))
            end
            return append!(lhs, rhs)

        @case (Ω::Number, ϕ) && if is_time_function(ϕ) end
            return [
                (t->exp(ϕ(t) * im)) => sum(put(ex.nsites, i=>matblock([0 Ω;0 0])) for i in 1:ex.nsites),
                (t->exp(-ϕ(t) * im)) => sum(put(ex.nsites, i=>matblock([0 0;Ω 0])) for i in 1:ex.nsites),
            ]

        @case (Ω::Vector{<:Number}, ϕ::Vector) && if is_time_function(ϕ) end
            lhs = map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
                (t->exp(ϕ_i(t) * im)) => put(ex.nsites, i=>matblock([0 Ω_i;0 0]))
            end

            rhs = map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
                (t->exp(-ϕ_i(t) * im)) => put(ex.nsites, i=>matblock([0 0;Ω_i 0]))
            end
            return append!(lhs, rhs)
        @case (Ω::Vector{<:Number}, ϕ) && if is_time_function(ϕ) end
            return [
                (t->exp(ϕ(t) * im)) => sum(put(ex.nsites, i=>matblock([0 Ω_i;0 0])) for (i, Ω_i) in enumerate(Ω)),
                (t->exp(-ϕ(t) * im)) => sum(put(ex.nsites, i=>matblock([0 0;Ω_i 0])) for (i, Ω_i) in enumerate(Ω)),
            ]

        # both time-dependent
        @case (Ω::Vector, ϕ::Vector) && if is_time_function(ϕ) && is_time_function(Ω) end
            lhs = map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
                (t->Ω_i(t) * exp(ϕ_i(t) * im)) => put(ex.nsites, i=>matblock([0 1;0 0]))
            end

            rhs = map(enumerate(zip(Ω, ϕ))) do (i, (Ω_i, ϕ_i))
                (t->Ω_i(t) * exp(-ϕ_i(t) * im)) => put(ex.nsites, i=>matblock([0 0;1 0]))
            end
            return append!(lhs, rhs)
        @case (Ω::Vector, ϕ) && if is_time_function(ϕ) && is_time_function(Ω) end
            lhs = map(enumerate(Ω)) do (i, Ω_i)
                (t->Ω_i(t) * exp(ϕ(t) * im)) => put(ex.nsites, i=>matblock([0 1;0 0]))
            end

            rhs = map(enumerate(Ω)) do (i, Ω_i)
                (t->Ω_i(t) * exp(-ϕ(t) * im)) => put(ex.nsites, i=>matblock([0 0;1 0]))
            end
            return append!(lhs, rhs)
        @case (Ω, ϕ::Vector) && if is_time_function(ϕ) && is_time_function(Ω) end
            lhs = map(enumerate(ϕ)) do (i, ϕ_i)
                (t->Ω(t) * exp(ϕ_i(t) * im)) => put(ex.nsites, i=>matblock([0 1;0 0]))
            end

            rhs = map(enumerate(ϕ)) do (i, ϕ_i)
                (t->Ω(t) * exp(-ϕ_i(t) * im)) => put(ex.nsites, i=>matblock([0 0;1 0]))
            end
            return append!(lhs, rhs)
        @case (Ω, ϕ) && if is_time_function(ϕ) && is_time_function(Ω) end
            return [
                (t->Ω(t) * exp(ϕ(t) * im)) => sum(put(ex.nsites, i=>matblock([0 1;0 0])) for i in 1:ex.nsites),
                (t->Ω(t) * exp(-ϕ(t) * im)) => sum(put(ex.nsites, i=>matblock([0 0;1 0])) for i in 1:ex.nsites),
            ]
        @case (Ω, ϕ)
            error("unexpected value for Ω and ϕ, got $(typeof((Ω, ϕ)))")
    end
end

"""
    Hamiltonian(::Type{Tv}, expr[, space=fullspace])

Create a `Hamiltonian` from hamiltonian expr that has matrix element of type `Tv`.
"""
function Hamiltonian(::Type{Tv}, ex::AbstractBlock, space::AbstractSpace=fullspace) where {Tv}
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
    push!(ts, mat(Tv, const_term, space))
    return Hamiltonian((fs..., ), (ts..., ))
end

function YaoBlocks.Optimise.to_basictypes(h::SumOfX)
    is_time_function(h.Ω) && error("cannot get matrix of a time-dependent operator")
    return if h.Ω isa Vector
        sum(h.Ω[i] * put(h.nsites, i=>X) for i in 1:h.nsites)
    else
        h.Ω * sum(put(h.nsites, i=>X) for i in 1:h.nsites)
    end
end

function YaoBlocks.Optimise.to_basictypes(h::Union{SumOfN, SumOfZ})
    is_time_function(h.Δ) && error("cannot get matrix of a time-dependent operator")
    op = h isa SumOfN ? ConstGate.P1 : Z
    return if h.Δ isa Vector
        sum(h.Δ[i] * put(h.nsites, i=>op) for i in 1:h.nsites)
    else
        h.Δ * sum(put(h.nsites, i=>op) for i in 1:h.nsites)
    end
end

function YaoBlocks.Optimise.to_basictypes(h::SumOfXPhase)
    (is_time_function(h.Ω) || is_time_function(h.ϕ)) &&
        error("cannot get matrix of a time-dependent operator")

    return @switch (h.Ω, h.ϕ) begin
        @case (::Vector, ::Vector)
            sum(h.Ω[i] * put(h.nsites, i=>XPhase(h.ϕ[i])) for i in 1:h.nsites)
        @case (Ω::Vector, ϕ)
            sum(Ω_i * put(h.nsites, i=>XPhase(ϕ)) for (i, Ω_i) in enumerate(Ω))
        @case (Ω, ϕ::Vector)
            sum(Ω * put(h.nsites, i=>XPhase(ϕ_i)) for (i, ϕ_i) in enumerate(ϕ))
        @case (Ω, ϕ)
            sum(Ω * put(h.nsites, i=>XPhase(ϕ)) for i in 1:h.nsites)
    end
end

function YaoBlocks.Optimise.to_basictypes(ex::RydInteract)
    nsites = length(ex.atoms)

    term = nothing
    for i in 1:nsites, j in 1:i-1
        x, y = ex.atoms[i], ex.atoms[j]
        h = ex.C / distance(x, y)^6 * kron(nsites, i=>ConstGate.P1, j=>ConstGate.P1)

        if isnothing(term)
            term = h
        else
            term += h
        end
    end
    return term
end

function emit_lowered(h)
    return YaoBlocks.Optimise.simplify(
        h; rules=[YaoBlocks.Optimise.to_basictypes]
    )
end
