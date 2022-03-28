YaoBlocks.print_block(io::IO, t::AbstractTerm) = print_expr(io, MIME"text/plain"(), t)
YaoBlocks.print_block(io::IO, t::XPhase) = print(io, "XPhase(", t.ϕ, ")")

Base.show(io::IO, ::Union{MIME"text/latex", MIME"application/x-latex"}, t::AbstractTerm) = print(io, latex_expr(t))
Base.show(io::IO, ::Union{MIME"text/latex", MIME"application/x-latex"}, t::Add) = print(io, latex_expr(t))

function latex_expr(t::Add)
    tex = latex_expr(t.list[1])

    for idx in 2:length(t.list)
        blk = t.list[idx]

        if blk isa Scale
            alpha = factor(blk)
            blk = content(blk)
            if alpha == -1
                tex = latexstring(tex, "-")
            else
                tex = latexstring(tex, "+", alpha, "\\cdot")
            end
        else
            tex = latexstring(tex, "+")
        end
        tex = latexstring(tex, latex_expr(blk))
    end
    return tex
end

function print_expr(io::IO, ::MIME"text/plain", t::RydInteract)
    C = round(t.C, sigdigits=3)
    print(io, "∑ $C/|r_i-r_j|^6 n_i n_j")
end

function latex_expr(t::RydInteract)
    C = round(t.C, sigdigits=3)
    return latexstring("\\sum \\frac{$C}{|r_i-r_j|^6} n_i n_j")
end

function print_expr(io::IO, ::MIME"text/plain", t::SumOfX)
    if t.Ω isa Number
        if isone(t.Ω)
            print(io, "∑ σ^x_i")
        else
            Ω = round(t.Ω, sigdigits=3)
            print(io, "$Ω ⋅ ∑ σ^x_i")
        end
    elseif t.Ω isa Vector
        print(io, "∑ Ω_i ⋅ σ^x")
    else
        print(io, "Ω(t) ⋅ ∑ σ^x_i")
    end
end

function latex_expr(t::SumOfX)
    if t.Ω isa Number
        if isone(t.Ω)
            tex = "\\sum σ^x_i"
        else
            Ω = round(t.Ω, sigdigits=3)
            tex = "$Ω ⋅ \\sum σ^x_i"
        end
    elseif t.Ω isa Vector
        tex = "\\sum Ω_i ⋅ σ^x"
    else
        tex = "Ω(t) ⋅ \\sum σ^x_i"
    end
    return latexstring(tex)
end

function print_expr(io::IO, ::MIME"text/plain", t::Union{SumOfN, SumOfZ})
    op = t isa SumOfN ? "n_i" : "σ^z_i"
    if t.Δ isa Number
        if isone(t.Δ)
            print(io, "∑ $op")
        else
            Δ = round(t.Δ, sigdigits=3)
            print(io, "$Δ ⋅ ∑ $op")
        end
    elseif t.Δ isa Vector
        print(io, "∑ Δ_i ⋅ $op")
    else
        print(io, "Δ(t) ⋅ ∑ $op")
    end
end

function latex_expr(t::Union{SumOfN, SumOfZ})
    op = t isa SumOfN ? "n_i" : "σ^z_i"
    if t.Δ isa Number
        if isone(t.Δ)
            tex = "\\sum $op"
        else
            Δ = round(t.Δ, sigdigits=3)
            tex = "$Δ ⋅ \\sum $op"
        end
    elseif t.Δ isa Vector
        tex = "\\sum Δ_i ⋅ $op"
    else
        tex = "Δ(t) ⋅ \\sum $op"
    end
    return latexstring(tex)
end

function print_expr(io::IO, ::MIME"text/plain", t::SumOfXPhase)
    @switch (t.Ω, t.ϕ) begin
        @case (Ω::Vector, ϕ::Vector)
            Ω = is_time_function(Ω) ? "Ω(t)_i" : "Ω_i"
            ϕ = is_time_function(ϕ) ? "ϕ(t)_i" : "ϕ_i"
            print(io, "∑ $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)")
        @case (Ω::Vector, ϕ)
            Ω = is_time_function(Ω) ? "Ω(t)_i" : "Ω_i"
            ϕ = is_time_function(ϕ) ? "ϕ(t)" : round(ϕ, sigdigits=3)
            print(io, "∑ $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)")
        @case (Ω, ϕ::Vector)
            Ω = is_time_function(Ω) ? "Ω(t)" : round(Ω, sigdigits=3)
            ϕ = is_time_function(ϕ) ? "ϕ(t)_i" : "ϕ_i"
            print(io, "$Ω ⋅ ∑ e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|")
        @case (Ω, ϕ)
            Ω = is_time_function(Ω) ? "Ω(t)" : round(Ω, sigdigits=3)
            ϕ = is_time_function(ϕ) ? "ϕ(t)" : round(ϕ, sigdigits=3)
            print(io, "$Ω ⋅ ∑ e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|")
    end
end

function latex_expr(t::SumOfXPhase)
    @switch (t.Ω, t.ϕ) begin
        @case (Ω::Vector, ϕ::Vector)
            Ω = is_time_function(Ω) ? "Ω(t)_i" : "Ω_i"
            ϕ = is_time_function(ϕ) ? "ϕ(t)_i" : "ϕ_i"
            tex = "\\sum $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)"
        @case (Ω::Vector, ϕ)
            Ω = is_time_function(Ω) ? "Ω(t)_i" : "Ω_i"
            ϕ = is_time_function(ϕ) ? "ϕ(t)" : round(ϕ, sigdigits=3)
            tex = "\\sum $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)"
        @case (Ω, ϕ::Vector)
            Ω = is_time_function(Ω) ? "Ω(t)" : round(Ω, sigdigits=3)
            ϕ = is_time_function(ϕ) ? "ϕ(t)_i" : "ϕ_i"
            tex = "$Ω ⋅ \\sum e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|"
        @case (Ω, ϕ)
            Ω = is_time_function(Ω) ? "Ω(t)" : round(Ω, sigdigits=3)
            ϕ = is_time_function(ϕ) ? "ϕ(t)" : round(ϕ, sigdigits=3)
            tex = "$Ω ⋅ \\sum e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|"
    end
    return tex
end

function Base.show(io::IO, ::MIME"text/plain", h::Hamiltonian)
    tab(n) = " "^(n + get(io, :indent, 0))
    println(io, "Hamiltonian")
    println(io, tab(2), "number of dynamic terms: ", count(x->x!==Base.one, h.fs))
    print(io, tab(2), "storage size: ")
    print(io, Base.format_bytes(sum(sizeof, h.ts)))
end
