distance(a, b) = sqrt(mapreduce(x->x^2, +, a .- b))

function YaoAPI.mat(::Type{T}, h::AbstractBlock, space::FullSpace) where T
    return YaoAPI.mat(T, h)
end

function YaoAPI.mat(::Type{T}, h::XPhase) where T
    return PermMatrix([2, 1], [exp(h.ϕ*im), exp(-h.ϕ*im)])
end

function YaoAPI.mat(::Type{T}, h::RydInteract) where T
    n = nqubits(h)
    values = Vector{T}(undef, 2^n)
    for idx in 1:2^n
        v = zero(T)
        for i in 1:n, j in 1:i-1
            if (readbit(idx, i) == 1) && (readbit(idx, j) == 1)
                r_i, r_j = h.atoms[i], h.atoms[j]
                alpha = h.C / distance(r_i, r_j)^6
                v += alpha
            end
        end
        values[idx] = v
    end
    return Diagonal(values)
end

function YaoAPI.mat(::Type{T}, h::SumOfX) where T
    is_time_function(h.Ω) && error("cannot get matrix of a time-dependent operator")
    ex = if h.Ω isa Vector
        sum(h.Ω[i] * put(h.nsites, i=>X) for i in 1:h.nsites)
    else
        h.Ω * sum(put(h.nsites, i=>X) for i in 1:h.nsites)
    end
    return mat(T, ex)
end

function YaoAPI.mat(::Type{T}, h::Union{SumOfN, SumOfZ}) where T
    is_time_function(h.Δ) && error("cannot get matrix of a time-dependent operator")
    op = h isa SumOfN ? ConstGate.P1 : Z
    ex = if h.Δ isa Vector
        sum(h.Δ[i] * put(h.nsites, i=>op) for i in 1:h.nsites)
    else
        h.Δ * sum(put(h.nsites, i=>op) for i in 1:h.nsites)
    end
    return mat(T, ex)
end

function YaoAPI.mat(::Type{T}, h::SumOfXPhase) where T
    (is_time_function(h.Ω) || is_time_function(h.ϕ)) &&
        error("cannot get matrix of a time-dependent operator")

    ex = @switch (h.Ω, h.ϕ) begin
        @case (::Vector, ::Vector)
            sum(h.Ω[i] * put(h.nsites, i=>XPhase(h.ϕ[i])) for i in 1:h.nsites)
        @case (Ω::Vector, ϕ)
            sum(Ω_i * put(h.nsites, i=>XPhase(ϕ)) for (i, Ω_i) in enumerate(Ω))
        @case (Ω, ϕ::Vector)
            sum(Ω * put(h.nsites, i=>XPhase(ϕ_i)) for (i, ϕ_i) in enumerate(ϕ))
        @case (Ω, ϕ)
            sum(Ω, put(h.nsites, i=>XPhase(ϕ)) for i in 1:h.nsites)
    end

    return mat(T, ex)
end
