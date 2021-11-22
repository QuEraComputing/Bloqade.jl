abstract type AbstractBravaisLattice end

struct BravaisLattice{D, T<:Real}
    bases::NTuple{D,NTuple{D,T}}
    sites::Vector{NTuple{D,T}}
end

function generate_lattice(bl::BravaisLattice{2,T}, m::Int, n::Int) where {T}
    locations = Tuple{T,T}[]  # we might want to avoid using `push!` later.
    for j=1:n
        for i=1:m
            baseloc = (i-1) .* bl.bases[1] .+ (j-1) .* bl.bases[2]
            for siteloc in bl.sites
                push!(locations, baseloc .+ siteloc)
            end
        end
    end
    return locations
end

bases = ((1.0, 0.0), (0.5, 0.5*sqrt(3)))
sites = [(0.0, 0.0), (0.5, 0.5/sqrt(3))]
honey = BravaisLattice(bases, sites)

using Viznet, Compose
import Cairo
locations = generate_lattice(honey, 5, 5)

function convert_locs(atoms::Vector{<:Tuple})
    xmin = minimum(x->x[1], atoms)
    ymin = minimum(x->x[2], atoms)
    xmax = maximum(x->x[1], atoms)
    ymax = maximum(x->x[2], atoms)
    scale = max(xmax-ymin, ymax-ymin)*1.1
    newlocs = map(x-> (x .- (xmin, ymin)) ./ scale .+ 0.05, atoms)
end

function viz_atoms(atoms::Vector{<:Tuple})
    line_style=compose(bondstyle(:default), stroke("black"))
    node_style=compose(context(), nodestyle(:default), stroke("black"), fill("white"), linewidth(0.5mm))
    text_style=textstyle(:default, fontsize(8pt))
    newlocs = convert_locs(atoms)
    canvas() do
        for (i, node) in enumerate(newlocs)
            node_style >> node
            text_style >> (node, "$i")
        end
    end
end
