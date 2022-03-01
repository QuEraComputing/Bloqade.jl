function get_matrix(::Type{Tv}, op::AbstractBlock, ::FullSpace) where Tv
    return mat(Tv, op)
end

function get_matrix(::Type{Tv}, op::AbstractBlock, space::Subspace) where Tv
    return mat(Tv, op, space)
end

function get_matrix(::Type{Tv}, op::AbstractTerm, space::AbstractSpace) where Tv
    return SparseMatrixCSC{Tv}(op, space)
end

struct ConstTermCache{FS <: Tuple, HS <: Tuple}
    fs::FS # time-dependent factors
    hs::HS # const terms
end

function storage_size(h::ConstTermCache)
    return sum(storage_size, h.hs)
end

# split const term and its dynamic prefactors from hamiltonian expr
function split_const_term(::Type{Tv}, h::Hamiltonian, space::AbstractSpace) where {Tv}
    fs, hs = [], []
    for t in h.terms, (f, h) in _split_term(Tv, t, space)
        push!(fs, f)
        # NOTE: we force converting blocks to a matrix as a workaround
        # of https://github.com/QuantumBFS/BQCESubroutine.jl/issues/37
        # so that we don't need to special case blocks to preallocate
        # the intermediate state for dstate.
        if h isa AbstractBlock
            error("unexpected block object")
        elseif h isa SparseMatrixCSC
            # always use CSR since it's faster in gemv
            push!(hs, transpose(SparseMatrixCSC(transpose(h))))
        else # other matrix type, e.g PermMatrix
            push!(hs, h)
        end
    end
    return ConstTermCache((fs...,), (hs...,))
end

function _split_term(::Type{Tv}, h::RydInteract, space::AbstractSpace) where {Tv}
    # TODO: actually implement it as Diagonal
    ((_const_param_, Diagonal(Vector(diag(SparseMatrixCSC{Tv}(h, space))))), )
end

function _split_term(::Type{Tv}, h::Negative, space::AbstractSpace) where {Tv}
    return map(_split_term(Tv, h.term, space)) do (f, h)
        f, -h
    end
end

_const_param_(t) = one(t)

function _split_term(::Type{Tv}, h::XTerm, space::AbstractSpace) where {Tv}
    n = nsites(h)
    @switch (h.Ωs, h.ϕs) begin
        @case (Ωs::ConstParamListType, ϕ::Number) || (Ωs::ConstParamListType, ::Nothing) || (Ω::Number, ϕ::Number) ||
        (Ω::Number, ::ConstParamListType) || (Ω::Number, ::Nothing)
            ((_const_param_, SparseMatrixCSC{Tv, Cint}(h, space)), )
        @case (Ωs::AbstractVector, ϕs::ConstParamListType) # directly apply is faster
            map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                x_phase = PermMatrix([2, 1], Tv[exp(ϕ * im), exp(-ϕ * im)])
                t->Ω(t)/2, get_matrix(Tv, put(n, i => matblock(x_phase)), space)
            end
        @case (Ωs::ConstParamListType, ϕs::ParamsList) # directly apply is faster
            op1 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω/2 * exp(ϕ(t) * im)), get_matrix(Tv, put(n, i => matblock(Tv[0 1;0 0])), space)
            end

            op2 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω/2 * exp(-ϕ(t) * im)), get_matrix(Tv, put(n, i => matblock(Tv[0 0;1 0])), space)
            end
            return (op1..., op2...)
        @case (Ωs::ParamsList, ϕs::ParamsList)
            op1 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω(t)/2 * exp(ϕ(t) * im)), get_matrix(Tv, put(n, i => matblock(Tv[0 1;0 0])), space)
            end

            op2 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω(t)/2 * exp(-ϕ(t) * im)), get_matrix(Tv, put(n, i => matblock(Tv[0 0;1 0])), space)
            end
            return (op1..., op2...)
        @case (Ωs::ConstParamListType, ϕ)
            op1 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω/2 * exp(ϕ(t) * im)), get_matrix(Tv, put(n, i => matblock(Tv[0 1;0 0])), space)
            end

            op2 = map(enumerate(zip(Ωs, ϕs))) do (i, (Ω, ϕ))
                t->(Ω/2 * exp(-ϕ(t) * im)), get_matrix(Tv, put(n, i => matblock(Tv[0 0;1 0])), space)
            end
            return (op1..., op2...)
        @case (Ωs::ParamsList, ::Nothing)
            map(enumerate(Ωs)) do (i, Ω)
                t->Ω(t)/2, get_matrix(Tv, put(n, i=>X), space)
            end
        @case (Ω::Number, ::ParamsList)
            op1 = map(enumerate(ϕs)) do (i, ϕ)
                t->(Ω/2 * exp(ϕ(t) * im)), get_matrix(Tv, put(n, i => matblock(Tv[0 1;0 0])), space)
            end

            op2 = map(enumerate(ϕs)) do (i, ϕ)
                t->(Ω/2 * exp(-ϕ(t) * im)), get_matrix(Tv, put(n, i => matblock(Tv[0 0;1 0])), space)
            end
            return (op1..., op2...)
        @case (Ω, ϕ::Number)
            A = get_matrix(Tv, sum(put(n, i=>matblock(Tv[0 1;0 0])) for i in 1:n), space)
            B = get_matrix(Tv, sum(put(n, i=>matblock(Tv[0 0;1 0])) for i in 1:n), space)
            return (t->Ω(t)/2 * exp(ϕ * im), A), (t->Ω(t)/2 * exp(-ϕ * im), B)
        @case (Ω::Number, ϕ)
            A = get_matrix(Tv, sum(put(n, i=>matblock(Tv[0 1;0 0])) for i in 1:n), space)
            B = get_matrix(Tv, sum(put(n, i=>matblock(Tv[0 0;1 0])) for i in 1:n), space)
            return (t->Ω/2 * exp(ϕ(t) * im), A), (t->Ω/2 * exp(-ϕ(t) * im), B)
        @case (Ω, ::Nothing) # no 1/2 in prefactor, it's in the matrix already
            return ((t->Ω(t), SparseMatrixCSC{Tv, Cint}(XTerm(n, 1.0), space)), )
        @case (Ω, ϕ)
            A = get_matrix(Tv, sum(put(n, i=>matblock(Tv[0 1;0 0])) for i in 1:n), space)
            B = get_matrix(Tv, sum(put(n, i=>matblock(Tv[0 0;1 0])) for i in 1:n), space)
            return (t->Ω(t)/2 * exp(ϕ(t) * im), A), (t->Ω(t)/2 * exp(-ϕ(t) * im), B)
    end
end

function _split_term(::Type{Tv}, h::NTerm, space::AbstractSpace) where {Tv}
    n = nsites(h)
    return if h.Δs isa ConstParamType
        M = Diagonal(Vector(diag(SparseMatrixCSC{Tv}(h, space))))
        ((_const_param_, M), )
    elseif h.Δs isa ParamsList
        return map(enumerate(h.Δs)) do (i, Δ)
            Δ, get_matrix(Tv, put(n, i=>Yao.ConstGate.P1), space)
        end
    else
        M = Diagonal(Vector(diag(SparseMatrixCSC{Tv}(NTerm(n, one(Tv)), space))))
        return ((h.Δs, M), )
    end
end
