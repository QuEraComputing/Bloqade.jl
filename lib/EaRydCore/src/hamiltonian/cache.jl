struct XPhase{T} <: PrimitiveBlock{2}
    theta::T
end

Yao.nqudits(::XPhase) = 1
Yao.mat(::Type{T}, g::XPhase) where T = PermMatrix([2, 1], T[exp(g.theta * im), exp(-g.theta * im)])

struct ConstTermCache{FS <: Tuple, HS <: Tuple}
    fs::FS # time-dependent factors
    hs::HS # const terms
end

function split_term(::Type{Tv}, h::RydbergHamiltonian, space::AbstractSpace) where {Tv}
    return (
        split_term(Tv, h.interact, space)...,
        split_term(Tv, h.xterm, space)...,
        split_term(Tv, h.nterm, space)...
    )
end

function split_term(::Type{Tv}, h::RydInteract, space::AbstractSpace) where {Tv}
    ((identity, SparseMatrixCSC{Tv}(h, space)), )
end

Yao.mat(::Type{T}, x::AbstractBlock, ::FullSpace) = Yao.mat(T, x)

function split_term(::Type{Tv}, h::XTerm, space::AbstractSpace) where {Tv}
    n = nsites(h)
    @switch (h.Ωs, h.ϕs) begin
        @case (Ωs::ConstParamListType, ϕ::Number) || (Ωs::ConstParamListType, ::Nothing) || (Ω::Number, ϕ::Number) ||
        (Ω::Number, ::ConstParamListType) || (Ω::Number, ::Nothing)
            ((identity, SparseMatrixCSC{Tv}(h, space)), )
        @case (Ωs::AbstractVector, ϕs::ConstParamListType) # directly apply is faster
            map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->Ω(t)/2, put(n, i => XPhase(ϕ))
            end
        @case (Ωs::ConstParamListType, ϕs::ParamsList) # directly apply is faster
            op1 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω/2 * exp(ϕ(t) * im)), put(n, i => matblock(Tv[0 1;0 0]))
            end

            op2 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω/2 * exp(-ϕ(t) * im)), put(n, i => matblock(Tv[0 0;1 0]))
            end
            return (op1..., op2...)
        @case (Ωs::ParamsList, ϕs::ParamsList)
            op1 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω(t)/2 * exp(ϕ(t) * im)), put(n, i => matblock(Tv[0 1;0 0]))
            end

            op2 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω(t)/2 * exp(-ϕ(t) * im)), put(n, i => matblock(Tv[0 0;1 0]))
            end
            return (op1..., op2...)
        @case (Ωs::ConstParamListType, ϕ)
            op1 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω/2 * exp(ϕ(t) * im)), put(n, i => matblock(Tv[0 1;0 0]))
            end

            op2 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω/2 * exp(-ϕ(t) * im)), put(n, i => matblock(Tv[0 0;1 0]))
            end
            return (op1..., op2...)
        @case (Ωs::ParamsList, ::Nothing)
            map(enumerate(Ωs)) do (i, Ω)
                t->Ω(t)/2, put(n, i=>X)
            end
        @case (Ω::Number, ::ParamsList)
            op1 = map(enumerate(ϕs)) do (i, ϕ)
                t->(Ω/2 * exp(ϕ(t) * im)), put(n, i => matblock(Tv[0 1;0 0]))
            end

            op2 = map(enumerate(ϕs)) do (i, ϕ)
                t->(Ω/2 * exp(-ϕ(t) * im)), put(n, i => matblock(Tv[0 0;1 0]))
            end
            return (op1..., op2...)
        @case (Ω, ϕ::Number)
            A = mat(Tv, sum(put(n, i=>matblock(Tv[0 1;0 0]))))
            B = mat(Tv, sum(put(n, i=>matblock(Tv[0 0;1 0]))))
            return (t->Ω(t)/2 * exp(ϕ * im), A), (t->Ω(t)/2 * exp(-ϕ * im), B)
        @case (Ω, ::Nothing)
            return ((t->Ω(t)/2, SparseMatrixCSC{Tv}(XTerm(n, 1.0), space)), )
        @case (Ω, ϕ)
            A = mat(Tv, sum(put(n, i=>matblock(Tv[0 1;0 0]))))
            B = mat(Tv, sum(put(n, i=>matblock(Tv[0 0;1 0]))))
            return (t->Ω(t)/2 * exp(ϕ(t) * im), A), (t->Ω(t)/2 * exp(-ϕ(t) * im), B)
    end
end

function split_term(::Type{Tv}, h::NTerm, space::AbstractSpace) where {Tv}
    n = nsites(h)
    return if h.Δs isa ConstParamType
        ((identity, -SparseMatrixCSC{Tv}(h, space)), )
    elseif h.Δs isa ParamsList
        return map(enumerate(h.Δs)) do (i, Δ)
            Δ, -put(n, i=>Yao.P1)
        end
    else
        return ((h.Δs, -SparseMatrixCSC{Tv}(NTerm(n, one(Tv)), space)), )
    end
end
