function latex_string(t::RydInteract)
    C = round(t.C / (2π); sigdigits=6)
    return "\\sum_{<ij>} \\frac{2π ⋅ $(C) MHz⋅μm}{|r_i - r_j|^6} n_i⋅n_j"
end

function latex_string(t::XTerm)
    n = nsites(t)
    @switch (t.Ωs, t.ϕs) begin
        @case (Ωs::AbstractVector{<:Number}, ϕs::AbstractVector{<:Number})
            tex = "\\sum_{k=1}^{$n} Ω_k (e^{ϕ_{k}i}|0⟩⟨1|_{k} + e^{-ϕ_{k}i}|1⟩⟨0|_{k})"
        @case (Ωs::AbstractVector, ϕs::AbstractVector{<:Number})
            tex = "\\sum_{k=1}^{$n} Ω_k(t) (e^{ϕ_{k}i}|0⟩⟨1|_{k} + e^{-ϕ_{k}i}|1⟩⟨0|_{k})"
        @case (Ωs::AbstractVector{<:Number}, ϕs::AbstractVector)
            tex = "\\sum_{k=1}^{$n} Ω_k (e^{ϕ_{k}(t)i}|0⟩⟨1|_{k} + e^{-ϕ_{k}(t)i}|1⟩⟨0|_{k})"
        @case (Ωs::AbstractVector, ϕs::AbstractVector)
            tex = "\\sum_{k=1}^{$n} Ω_k(t) (e^{ϕ_{k}(t)i}|0⟩⟨1|_{k} + e^{-ϕ_{k}(t)i}|1⟩⟨0|_{k})"
        @case (Ωs::AbstractVector{<:Number}, ϕ::Number)
            ϕ = round(ϕ;sigdigits=4)
            tex = "\\sum_{k=1}^{$n} Ω_k (e^{$(ϕ)i}|0⟩⟨1|_{k} + e^{-$(ϕ)i}|1⟩⟨0|_{k})"
        @case (Ωs::AbstractVector{<:Number}, ::Nothing)
            tex = "\\sum_{k=1}^{$n} Ω_k σ^x_k"
        @case (Ωs::AbstractVector{<:Number}, ϕ)
            tex = "\\sum_{k=1}^{$n} Ω_k (e^{ϕ(t)i}|0⟩⟨1|_{k} + e^{-ϕ(t)i}|1⟩⟨0|_{k})"
        @case (Ωs::AbstractVector, ::Nothing)
            tex = "\\sum_{k=1}^{$n} Ω_k(t) σ^x_k"
        @case (Ω::Number, ϕ::Number)
            Ω, ϕ = round(Ω;sigdigits=6), round(ϕ;sigdigits=4)
            tex = "$(Ω)\\sum_{k=1}^{$n} (e^{$(ϕ)i}|0⟩⟨1|_{k} + e^{-$(ϕ)i}|1⟩⟨0|_{k})"
        @case (Ω::Number, ::AbstractVector{<:Number})
            tex = "$(Ω)\\sum_{k=1}^{$n} (e^{ϕ_{k}i}|0⟩⟨1|_{k} + e^{-ϕ_{k}i}|1⟩⟨0|_{k})"
        @case (Ω::Number, ::AbstractVector)
            tex = "$(Ω)\\sum_{k=1}^{$n} (e^{ϕ_{k}(t)i}|0⟩⟨1|_{k} + e^{-ϕ_{k}(t)i}|1⟩⟨0|_{k})"
        @case (Ω::Number, ::Nothing)
            tex = "$Ω\\sum_{k=1}^{$n} σ^x_k"
        @case (Ω, ::Nothing)
            tex = "Ω(t)\\sum_{k=1}^{$n} σ^x_k"
        @case (Ω, ϕ)
            tex = "Ω(t)\\sum_{k=1}^{$n} (e^{ϕ(t)i}|0⟩⟨1|_{k} + e^{-ϕ(t)i}|1⟩⟨0|_{k})"
    end
    return tex
end

function latex_string(t::NTerm)
    n = nsites(t)
    return if t.Δs isa AbstractVector{<:Number}
        "\\sum_{k=1}^{$n} Δ_k ⋅ n_k"
    elseif t.Δs isa AbstractVector
        "\\sum_{k=1}^{$n} Δ(t)_k ⋅ n_k"
    elseif t.Δs isa Number
        "$(t.Δs)\\sum_{k=1}^{$n} n_k"
    else
        "Δ(t)\\sum_{k=1}^{$n} n_k"
    end
end

function latex_string(t::Negative)
    return "-" * latex_string(t.term)
end

function latex_string(t::Hamiltonian)
    s = latex_string(first(t.terms))
    for term in t.terms[2:end]
        if term isa Negative
            s *= " " * latex_string(term)
        else
            s *= " + " * latex_string(term)
        end
    end
    return s
end

function Base.show(io::IO, ::MIME"text/latex", t::AbstractTerm)
    print(io, LaTeXStrings.latexstring(latex_string(t)))
end

# function Base.show(io::IO, ::MIME"text/latex", t::Hamiltonian)
# end

# function Base.show(io::IO, ::MIME"text/plain", t::AbstractTerm)
#     indent = get(io, :indent, 0)
#     println(io, " "^indent, nameof(typeof(t)))
#     print_term(IOContext(io, :indent=>indent + 1), t)
# end
