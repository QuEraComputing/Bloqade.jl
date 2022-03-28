module Op

using YaoBlocks
const n = YaoBlocks.ConstGate.P1

end

function rydberg_h(atom_positions; C::Real=2π * 862690, Ω=nothing, ϕ=nothing, Δ=nothing)
    positions = map(atom_positions) do pos
        (pos..., )
    end

    nsites = length(positions)
    term = RydInteract(positions, C)

    Ω = div_by_two(Ω)

    if !isnothing(Ω) && !isnothing(ϕ)
        term += SumOfXPhase(nsites, Ω, ϕ)
    elseif !isnothing(Ω) && isnothing(ϕ)
        term += SumOfX(nsites, Ω)
    end

    if !isnothing(Δ)
        term -= SumOfN(nsites, Δ)
    end

    return term
end

function div_by_two(Ω)
    isnothing(Ω) && return
    if !is_time_function(Ω)
        return Ω ./ 2
    end

    return if Ω isa Vector
        map(Ω) do Ω_i
            t->Ω_i(t) / 2
        end
    else
        t->Ω(t)/2
    end
end

attime(t::Real) = h->attime(h, t)

function attime(h::AbstractBlock, t::Real)
    blks = map(subblocks(h)) do blk
        attime(blk, t)
    end
    return chsubblocks(h, blks)
end

# Yao Blocks cannot take time-dependent function
attime(h::PrimitiveBlock, ::Real) = h
function attime(h::SumOfX, t::Real)
    if h.Ω isa Vector
        SumOfX(h.nsites, map(x->x(t), h.Ω))
    else
        SumOfX(h.nsites, h.Ω)
    end
end

function attime(h::SumOfXPhase, t::Real)
    if h.Ω isa Vector
        Ω = map(x->x(t), h.Ω)
    else
        Ω = h.Ω(t)
    end

    if h.ϕ isa Vector
        ϕ = map(x->x(t), h.ϕ)
    else
        ϕ = h.ϕ(t)
    end

    return SumOfXPhase(h.nsites, Ω, ϕ)
end

function attime(h::Union{SumOfZ, SumOfN}, t::Real)
    if h.Δ isa Vector
        Δ = map(x->x(t), h.Δ)
    else
        Δ = h.Δ(t)
    end
    return typeof(h)(h.nsites, Δ)
end
