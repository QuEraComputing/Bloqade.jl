function lower_expr(ex::Add)
    list = []
    for t in subblocks(ex)
        append!(list, lower_expr(t))
    end
    return list
end

function lower_expr(ex::Scale)
    return map(lower_expr(content(ex))) do (f, t)
        f => factor(ex) * t
    end
end

function lower_expr(ex::SumOfX)
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

function lower_expr(ex::Union{SumOfN, SumOfZ})
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

function lower_expr(ex::SumOfXPhase)
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


function Hamiltonian(::Type{Tv}, ex::AbstractBlock, space::AbstractSpace=fullspace) where {Tv}
    fs, ts = [], []
    for (f, op) in lower_expr(ex)
        push!(fs, f)
        push!(ts, mat(Tv, t, space))
    end
    return Hamiltonian((fs..., ), (ts..., ))
end
