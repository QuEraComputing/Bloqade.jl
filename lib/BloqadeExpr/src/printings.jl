YaoBlocks.print_block(io::IO, t::AbstractTerm) = print_expr(io, MIME"text/plain"(), t)
YaoBlocks.print_block(io::IO, t::RydbergHamiltonian) = print_expr(io, MIME"text/plain"(), t)
YaoBlocks.print_block(io::IO, t::XPhase) = print(io, "XPhase(", t.ϕ, ")")

Base.show(io::IO, ::Union{MIME"text/latex",MIME"application/x-latex"}, t::AbstractTerm) =
    print(io, latexstring(latex_expr(t)))
Base.show(io::IO, ::Union{MIME"text/latex",MIME"application/x-latex"}, t::RydbergHamiltonian) =
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

function print_expr(io::IO, ::MIME"text/plain", t::RydInteract{D}) where D
    C = t.C / 2π
    n = ceil(log10(C))-1
    C = round(C / 10^n, digits = 3)
    str_op = (D == 2 ? "n_i n_j" : "n^r_i n^r_j")
    return print(io, "∑ 2π ⋅ $(C)e$n/|x_i-x_j|^6 ", str_op)
end

function latex_expr(t::RydInteract{D}) where D
    C = t.C / 2π
    n = ceil(log10(C))
    C = round(C / 10^n, digits = 3)
    str_op = (D == 2 ? "n_i n_j" : "n^r_i n^r_j")
    return "\\sum \\frac{2π \\cdot $C\\times 10^{$n}}{|x_i-x_j|^6} " * str_op
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

function print_expr(io::IO, ::MIME"text/plain", t::SumOfXTypes)
    op = (t isa SumOfX ? "σ^x_i" : (t isa SumOfX_01 ? "σ^{x,hf}_i" : "σ^{x,r}_i"))
    if t.Ω isa Number
        print(io, pretty_number(t.Ω), "∑ ", op)
    elseif t.Ω isa Vector
        print(io, "∑ Ω_i ⋅ ", op)
    else
        print(io, "Ω(t) ⋅ ∑ ", op)
    end
end

function latex_expr(t::SumOfXTypes)
    op = (t isa SumOfX ? "σ^x_i" : (t isa SumOfX_01 ? "σ^{x,hr}_i" : "σ^{x,r}_i"))
    if t.Ω isa Number
        tex = pretty_number(t.Ω) * "\\sum " * op
    elseif t.Ω isa Vector
        tex = "\\sum Ω_i ⋅ " * op
    else
        tex = "Ω(t) ⋅ \\sum " * op
    end
    return tex
end

function print_expr(io::IO, ::MIME"text/plain", t::Union{SumOfN, SumOfN_1, SumOfN_r, SumOfZ, SumOfZ_01, SumOfZ_1r})
    op = if t isa SumOfN
        "n_i"
    elseif t isa SumOfN_1 
        "n^{hf}_i"
    elseif t isa SumOfN_r 
        "n^r_i"
    elseif t isa SumOfZ
        "σ^z_i"
    elseif t isa SumOfZ_01
        "σ^{z,hf}_i"
    elseif t isa SumOfZ_1r
        "σ^{z,r}_i"
    end
    if t.Δ isa Number
        print(io, pretty_number(t.Δ), "∑ $op")
    elseif t.Δ isa Vector
        print(io, "∑ Δ_i ⋅ $op")
    else
        print(io, "Δ(t) ⋅ ∑ $op")
    end
end

function latex_expr(t::Union{SumOfN, SumOfN_1, SumOfN_r, SumOfZ, SumOfZ_01, SumOfZ_1r}) 
    op = if t isa SumOfN
        "n_i"
    elseif t isa SumOfN_1 
        "n^{hf}_i"
    elseif t isa SumOfN_r 
        "n^r_i"
    elseif t isa SumOfZ
        "σ^z_i"
    elseif t isa SumOfZ_01
        "σ^{z,hf}_i"
    elseif t isa SumOfZ_1r
        "σ^{z,r}_i"
    end
    if t.Δ isa Number
        tex = pretty_number(t.Δ) * "\\sum $op"
    elseif t.Δ isa Vector
        tex = "\\sum Δ_i ⋅ $op"
    else
        tex = "Δ(t) ⋅ \\sum $op"
    end
    return tex
end

function print_expr(io::IO, ::MIME"text/plain", t::SumOfXPhaseTypes)
    if t isa SumOfXPhase
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
    elseif t isa SumOfXPhase_01
        @switch (t.Ω, t.ϕ) begin
            @case (Ω::Vector, ϕ::Vector)
            Ω = is_const_param(Ω) ? "Ω^{hf}_i" : "Ω^{hf}(t)_i"
            ϕ = is_const_param(ϕ) ? "ϕ^{hf}_i" : "ϕ^{hf}(t)_i"
            print(io, "∑ $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)")
            @case (Ω::Vector, ϕ)
            Ω = is_const_param(Ω) ? "Ω^{hf}_i" : "Ω^{hf}(t)_i"
            ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ^{hf}(t)"
            print(io, "∑ $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)")
            @case (Ω, ϕ::Vector)
            Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω^{hf}(t) ⋅"
            ϕ = is_const_param(ϕ) ? "ϕ^{hf}_i" : "ϕ^{hf}(t)_i"
            print(io, Ω, "∑ e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|")
            @case (Ω, ϕ)
            Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω^{hf}(t) ⋅"
            ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ^{hf}(t)"
            print(io, Ω, "∑ e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|")
        end
    elseif t isa SumOfXPhase_1r
        @switch (t.Ω, t.ϕ) begin
            @case (Ω::Vector, ϕ::Vector)
            Ω = is_const_param(Ω) ? "Ω^r_i" : "Ω^r(t)_i"
            ϕ = is_const_param(ϕ) ? "ϕ^r_i" : "ϕ^r(t)_i"
            print(io, "∑ $Ω ⋅ (e^{$ϕ ⋅ im} |1⟩⟨r| + e^{-$ϕ ⋅ im} |r⟩⟨1|)")
            @case (Ω::Vector, ϕ)
            Ω = is_const_param(Ω) ? "Ω^r_i" : "Ω^r(t)_i"
            ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ^r(t)"
            print(io, "∑ $Ω ⋅ (e^{$ϕ ⋅ im} |1⟩⟨r| + e^{-$ϕ ⋅ im} |r⟩⟨1|)")
            @case (Ω, ϕ::Vector)
            Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω^r(t) ⋅"
            ϕ = is_const_param(ϕ) ? "ϕ^r_i" : "ϕ^r(t)_i"
            print(io, Ω, "∑ e^{$ϕ ⋅ im} |1⟩⟨r| + e^{-$ϕ ⋅ im} |r⟩⟨1|")
            @case (Ω, ϕ)
            Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω^r(t) ⋅"
            ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ^r(t)"
            print(io, Ω, "∑ e^{$ϕ ⋅ im} |1⟩⟨r| + e^{-$ϕ ⋅ im} |r⟩⟨1|")
        end
    end
end

function latex_expr(t::SumOfXPhaseTypes)
    if t isa SumOfXPhase
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
    elseif t isa SumOfXPhase_01
        @switch (t.Ω, t.ϕ) begin
            @case (Ω::Vector, ϕ::Vector)
            Ω = is_const_param(Ω) ? "Ω^{hf}_i" : "Ω^{hf}(t)_i"
            ϕ = is_const_param(ϕ) ? "ϕ^{hf}_i" : "ϕ^{hf}(t)_i"
            tex = "\\sum $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)"
            @case (Ω::Vector, ϕ)
            Ω = is_const_param(Ω) ? "Ω^{hf}_i" : "Ω^{hf}(t)_i"
            ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ^{hf}(t)"
            tex = "\\sum $Ω ⋅ (e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|)"
            @case (Ω, ϕ::Vector)
            Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω^{hf}(t) ⋅ "
            ϕ = is_const_param(ϕ) ? "ϕ^{hf}_i" : "ϕ^{hf}(t)_i"
            tex = Ω * "\\sum e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|"
            @case (Ω, ϕ)
            Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω^{hf}(t) ⋅ "
            ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ^{hf}(t)"
            tex = Ω * "\\sum e^{$ϕ ⋅ im} |0⟩⟨1| + e^{-$ϕ ⋅ im} |1⟩⟨0|"
        end
    elseif t isa SumOfXPhase_1r
        @switch (t.Ω, t.ϕ) begin
            @case (Ω::Vector, ϕ::Vector)
            Ω = is_const_param(Ω) ? "Ω^r_i" : "Ω^r(t)_i"
            ϕ = is_const_param(ϕ) ? "ϕ^r_i" : "ϕ^r(t)_i"
            tex = "\\sum $Ω ⋅ (e^{$ϕ ⋅ im} |1⟩⟨r| + e^{-$ϕ ⋅ im} |r⟩⟨1|)"
            @case (Ω::Vector, ϕ)
            Ω = is_const_param(Ω) ? "Ω^r_i" : "Ω^r(t)_i"
            ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ^r(t)"
            tex = "\\sum $Ω ⋅ (e^{$ϕ ⋅ im} |1⟩⟨r| + e^{-$ϕ ⋅ im} |r⟩⟨1|)"
            @case (Ω, ϕ::Vector)
            Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω^r(t) ⋅ "
            ϕ = is_const_param(ϕ) ? "ϕ^r_i" : "ϕ^r(t)_i"
            tex = Ω * "\\sum e^{$ϕ ⋅ im} |1⟩⟨r| + e^{-$ϕ ⋅ im} |r⟩⟨1|"
            @case (Ω, ϕ)
            Ω = is_const_param(Ω) ? pretty_number(Ω) : "Ω^r(t) ⋅ "
            ϕ = is_const_param(ϕ) ? round(ϕ, sigdigits = 3) : "ϕ^r(t)"
            tex = Ω * "\\sum e^{$ϕ ⋅ im} |1⟩⟨r| + e^{-$ϕ ⋅ im} |r⟩⟨1|"
        end
    end
    return tex
end

function latex_expr(h::Union{RydbergHamiltonian, RydbergHamiltonian3})
    return latex_expr(add_terms(h))
end

function print_expr(io::IO, ::MIME"text/plain", h::Union{RydbergHamiltonian, RydbergHamiltonian3})
    print(io,add_terms(h))
end

function Base.show(io::IO, ::MIME"text/plain", h::Hamiltonian)
    tab(n) = " "^(n + get(io, :indent, 0))
    println(io, tab(0), "Hamiltonian")
    println(io, tab(2), "number of dynamic terms: ", count(x -> x !== Base.one, h.fs))
    print(io, tab(2), "storage size: ")
    return print(io, Base.format_bytes(sum(sizeof, h.ts)))
end
