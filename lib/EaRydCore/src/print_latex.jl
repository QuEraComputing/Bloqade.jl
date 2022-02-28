function latex_string(t::RydInteract)
    C = round(t.C / (2π); sigdigits=6)
    return "\\sum_{<ij>} \\frac{2π ⋅ $(C) MHz⋅μm^6}{|r_i - r_j|^6} n_i⋅n_j"
end

const ParamsList{T} = Union{AbstractVector{T}, NTuple{N, T} where N}

function latex_string(t::XTerm)
    n = nsites(t)
    @switch (t.Ωs, t.ϕs) begin
        @case (Ωs::ConstParamListType, ϕs::ConstParamListType)
            tex = "\\sum_{k=1}^{$n} \\frac{Ω_k}{2} (e^{ϕ_{k}i}|0⟩⟨1|_{k} + e^{-ϕ_{k}i}|1⟩⟨0|_{k})"
        @case (Ωs::AbstractVector, ϕs::ConstParamListType)
            tex = "\\sum_{k=1}^{$n} \\frac{Ω_k(t)}{2} (e^{ϕ_{k}i}|0⟩⟨1|_{k} + e^{-ϕ_{k}i}|1⟩⟨0|_{k})"
        @case (Ωs::ConstParamListType, ϕs::ParamsList)
            tex = "\\sum_{k=1}^{$n} \\frac{Ω_k}{2} (e^{ϕ_{k}(t)i}|0⟩⟨1|_{k} + e^{-ϕ_{k}(t)i}|1⟩⟨0|_{k})"
        @case (Ωs::ParamsList, ϕs::ParamsList)
            tex = "\\sum_{k=1}^{$n} \\frac{Ω_k(t)}{2} (e^{ϕ_{k}(t)i}|0⟩⟨1|_{k} + e^{-ϕ_{k}(t)i}|1⟩⟨0|_{k})"
        @case (Ωs::ConstParamListType, ϕ::Number)
            ϕ = round(ϕ;sigdigits=4)
            tex = "\\sum_{k=1}^{$n} \\frac{Ω_k}{2} (e^{$(ϕ)i}|0⟩⟨1|_{k} + e^{-$(ϕ)i}|1⟩⟨0|_{k})"
        @case (Ωs::ConstParamListType, ::Nothing)
            tex = "\\sum_{k=1}^{$n} \\frac{Ω_k}{2} σ^x_k"
        @case (Ωs::ConstParamListType, ϕ)
            tex = "\\sum_{k=1}^{$n} \\frac{Ω_k}{2} (e^{ϕ(t)i}|0⟩⟨1|_{k} + e^{-ϕ(t)i}|1⟩⟨0|_{k})"
        @case (Ωs::ParamsList, ::Nothing)
            tex = "\\sum_{k=1}^{$n} \\frac{Ω_k(t)}{2} σ^x_k"
        @case (Ω::Number, ϕ::Number)
            Ω, ϕ = round(Ω;sigdigits=6), round(ϕ;sigdigits=4)
            tex = "\\frac{$(Ω)}{2}\\sum_{k=1}^{$n} (e^{$(ϕ)i}|0⟩⟨1|_{k} + e^{-$(ϕ)i}|1⟩⟨0|_{k})"
        @case (Ω::Number, ::ConstParamListType)
            Ω = round(Ω;sigdigits=6)
            tex = "\\frac{$(Ω)}{2}\\sum_{k=1}^{$n} (e^{ϕ_{k}i}|0⟩⟨1|_{k} + e^{-ϕ_{k}i}|1⟩⟨0|_{k})"
        @case (Ω::Number, ::ParamsList)
            Ω = round(Ω;sigdigits=6)
            tex = "\\frac{$(Ω)}{2}\\sum_{k=1}^{$n} (e^{ϕ_{k}(t)i}|0⟩⟨1|_{k} + e^{-ϕ_{k}(t)i}|1⟩⟨0|_{k})"
        @case (Ω::Number, ::Nothing)
            Ω = round(Ω;sigdigits=6)
            tex = "\\frac{$(Ω)}{2}\\sum_{k=1}^{$n} σ^x_k"
        @case (Ω, ϕ::Number)
            ϕ = round(ϕ;sigdigits=4)
            tex = "\\frac{Ω(t)}{2}\\sum_{k=1}^{$n} (e^{$(ϕ)i}|0⟩⟨1|_{k} + e^{-$(ϕ)i}|1⟩⟨0|_{k})"
        @case (Ω, ::Nothing)
            tex = "\\frac{Ω(t)}{2}\\sum_{k=1}^{$n} σ^x_k"
        @case (Ω, ϕ)
            tex = "\\frac{Ω(t)}{2}\\sum_{k=1}^{$n} (e^{ϕ(t)i}|0⟩⟨1|_{k} + e^{-ϕ(t)i}|1⟩⟨0|_{k})"
    end
    return tex
end

function latex_string(t::NTerm)
    n = nsites(t)
    return if t.Δs isa ConstParamListType
        "\\sum_{k=1}^{$n} Δ_k ⋅ n_k"
    elseif t.Δs isa ParamsList
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
