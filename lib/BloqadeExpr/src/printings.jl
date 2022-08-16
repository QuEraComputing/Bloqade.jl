YaoBlocks.print_block(io::IO, t::AbstractTerm) = print_expr(io, MIME"text/plain"(), t)
YaoBlocks.print_block(io::IO, t::XPhase) = print(io, "XPhase(", t.ϕ, ")")

Base.show(io::IO, ::Union{MIME"text/latex",MIME"application/x-latex"}, t::AbstractTerm) =
    print(io, latexstring(latex_expr(t)))
Base.show(io::IO, ::Union{MIME"text/latex",MIME"application/x-latex"}, t::Add) = print(io, latexstring(latex_expr(t)))

function latex_expr(t::Add)
    tex = latex_expr(t.list[1])

    for idx in 2:length(t.list)
        blk = t.list[idx]

        if blk isa Scale
            alpha = factor(blk)
            blk = content(blk)
            if alpha == -1
                tex *= "-"
            else
                tex = string(tex, "+", alpha, "\\cdot")
            end
        else
            tex = string(tex, "+")
        end
        tex = string(tex, latex_expr(blk))
    end
    return tex
end

function latex_expr(t::YaoBlocks.Scale)
    if isone(factor(t))
        latex_expr(content(t))
    else
        string(factor(t), "\\cdot", latex_expr(content(t)))
    end
end

function print_expr(io::IO, ::MIME"text/plain", t::RydInteract)
    C = t.C / 2π
    n = ceil(log10(C))
    C = round(C / 10^(n - 1), digits = 3)
    return print(io, "∑ 2π ⋅ $(C)e$n/|r_i-r_j|^6 n_i n_j")
end

function latex_expr(t::RydInteract)
    C = t.C / 2π
    n = ceil(log10(C))
    C = round(C / 10^n, digits = 3)
    return "\\sum \\frac{2π \\cdot $C\\times 10^{$n}}{|r_i-r_j|^6} n_i n_j"
end

function pretty_number(x)
    return if isone(x)
        ""
    elseif isone(x / 2π)
        "2π ⋅ "
    else
        x = round(x / 2π, sigdigits = 3)
        "2π ⋅ $x ⋅ "
    end
end

function print_expr(io::IO, ::MIME"text/plain", t::SumOfX)
    if t.Ω isa Number
        print(io, pretty_number(t.Ω), "∑ σ^x_i")
    elseif t.Ω isa Vector
        print(io, "∑ Ω_i ⋅ σ^x")
    else
        print(io, "Ω(t) ⋅ ∑ σ^x_i")
    end
end

function latex_expr(t::SumOfX)
    if t.Ω isa Number
        tex = pretty_number(t.Ω) * "\\sum σ^x_i"
    elseif t.Ω isa Vector
        tex = "\\sum Ω_i ⋅ σ^x"
    else
        tex = "Ω(t) ⋅ \\sum σ^x_i"
    end
    return tex
end

function print_expr(io::IO, ::MIME"text/plain", t::Union{SumOfN,SumOfZ})
    op = t isa SumOfN ? "n_i" : "σ^z_i"
    if t.Δ isa Number
        print(io, pretty_number(t.Δ), "∑ $op")
    elseif t.Δ isa Vector
        print(io, "∑ Δ_i ⋅ $op")
    else
        print(io, "Δ(t) ⋅ ∑ $op")
    end
end

function latex_expr(t::Union{SumOfN,SumOfZ})
    op = t isa SumOfN ? "n_i" : "σ^z_i"
    if t.Δ isa Number
        tex = pretty_number(t.Δ) * "\\sum $op"
    elseif t.Δ isa Vector
        tex = "\\sum Δ_i ⋅ $op"
    else
        tex = "Δ(t) ⋅ \\sum $op"
    end
    return tex
end

function print_expr(io::IO, ::MIME"text/plain", t::SumOfXPhase)
    @switch (t.Ω, t.ϕ) begin
        @case (Ω::Vector, ϕ::Vector)
        Ω = is_const_param(Ω) ? "Ω_i" : "Ω(t)_i"
        ϕ = is_const_param(ϕ) ? "ϕ_i" : "ϕ(t)_i"
        print(io, "∑ $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)")
        @case (Ω::Vector, ϕ)
        Ω = is_const_param(Ω) ? "Ω_i" : "Ω(t)_i"
        ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ(t)"
        print(io, "∑ $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)")
        @case (Ω, ϕ::Vector)
        Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω(t) ⋅"
        ϕ = is_const_param(ϕ) ? "ϕ_i" : "ϕ(t)_i"
        print(io, Ω, "∑ e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|")
        @case (Ω, ϕ)
        Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω(t) ⋅"
        ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ(t)"
        print(io, Ω, "∑ e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|")
    end
end

function latex_expr(t::SumOfXPhase)
    @switch (t.Ω, t.ϕ) begin
        @case (Ω::Vector, ϕ::Vector)
        Ω = is_const_param(Ω) ? "Ω_i" : "Ω(t)_i"
        ϕ = is_const_param(ϕ) ? "ϕ_i" : "ϕ(t)_i"
        tex = "\\sum $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)"
        @case (Ω::Vector, ϕ)
        Ω = is_const_param(Ω) ? "Ω_i" : "Ω(t)_i"
        ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ(t)"
        tex = "\\sum $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)"
        @case (Ω, ϕ::Vector)
        Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω(t) ⋅ "
        ϕ = is_const_param(ϕ) ? "ϕ_i" : "ϕ(t)_i"
        tex = Ω * "\\sum e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|"
        @case (Ω, ϕ)
        Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω(t) ⋅ "
        ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ(t)"
        tex = Ω * "\\sum e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|"
    end
    return tex
end

function print_expr(io::IO, ::MIME"text/plain", h::RydbergHamiltonian)
    term = h.RydbergTerm

    if h.RabiTerm.Ω != 0
        term += h.RabiTerm
    end

    if h.DetuningTerm.Δ != 0
        term -= h.DetuningTerm
    end

    print(io,latex_expr(term))
end

function Base.show(io::IO, ::MIME"text/plain", h::Hamiltonian)
    tab(n) = " "^(n + get(io, :indent, 0))
    println(io, tab(0), "Hamiltonian")
    println(io, tab(2), "number of dynamic terms: ", count(x -> x !== Base.one, h.fs))
    print(io, tab(2), "storage size: ")
    return print(io, Base.format_bytes(sum(sizeof, h.ts)))
end
