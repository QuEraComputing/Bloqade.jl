struct Rescaler{T}
    xmin::T
    xmax::T
    ymin::T
    ymax::T
    pad::T
end

getscale(r::Rescaler) = 1/(r.xmax-r.xmin+r.pad)

function (r::Rescaler{T})(x; dims=(1,2)) where T
    xmin, ymin, xmax, ymax, pad = r.xmin, r.ymin, r.xmax, r.ymax, r.pad
    scale = getscale(r)
    if dims == (1,2)
        return (x[1]-xmin+0.5*pad, ymax+0.5*pad-x[2]) .* scale
    elseif dims == 1
        return (x - xmin + 0.5*pad) * scale
    elseif dims == 2
        return (ymax + 0.5*pad - x) * scale
    else
        throw(ArgumentError("dims should be (1,2), 1 or 2."))
    end
end

function get_rescaler(atoms::AbstractVector{<:Tuple}, pad)
    xmin = minimum(x->x[1], atoms)
    ymin = minimum(x->x[2], atoms)
    xmax = maximum(x->x[1], atoms)
    ymax = maximum(x->x[2], atoms)
    return Rescaler(xmin, xmax, ymin, ymax, pad)
end

default_node_style(scale; color="white") = compose(context(), Viznet.nodestyle(:default, r=0.15cm*scale), Compose.stroke("black"), fill(color), linewidth(0.3mm*scale))
default_text_style(scale) = Viznet.textstyle(:default, fontsize(4pt*scale))
default_bond_style(scale; color="black") = Viznet.bondstyle(:default, Compose.stroke(color), linewidth(0.3mm*scale))
default_line_style_grid(scale) = Viznet.bondstyle(:default, Compose.stroke("#AAAAAA"), linewidth(0.3mm*scale); dashed=true)
default_blockade_style(scale, blockade_radius) = Viznet.nodestyle(:circle, Compose.fill("transparent"), r=blockade_radius, Compose.stroke("black"), strokedash([0.5mm*scale, 0.5mm*scale]), linewidth(0.3mm*scale))

"""
    viz_atoms([io::Union{IO,AbstractString}, ]atoms::AtomList;
        colors=fill("white", count(maskedgrid.mask)), blockade_radius=0,
        texts=["1", "2", ...],
        blockade_radius=0;
        format=PNG,
        blockade_style="none",
        bond_color="black",
        )

Plots `atoms` with colors specified by `colors` and texts specified by `texts`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.
If you want to write this image to the disk without displaying it in a frontend, please try

```julia
julia> open("test.png", "w") do f
            viz_atoms(f, generate_sites(SquareLattice(), 5, 5))
       end
```

The `format` keyword argument can also be `Compose.SVG` or `Compose.PDF`.
"""
function viz_atoms(io, atoms::AtomList; format=PNG,
        colors=fill("white", length(atoms)),
        blockade_radius=0,
        texts = ["$i" for i=1:length(atoms)],
        kwargs...)
    img, (dx, dy) = img_atoms(atoms; colors=colors, blockade_radius=blockade_radius, texts=texts, config=LatticeDisplayConfig(; kwargs...))
    Compose.draw(format(io, dx, dy), img)
    return
end
function viz_atoms(atoms::AtomList;
        colors=fill("white", length(atoms)),
        blockade_radius=0,
        texts = ["$i" for i=1:length(atoms)],
        kwargs...)
    img, (dx, dy) = img_atoms(atoms; colors=colors, blockade_radius=blockade_radius, texts=texts, config=LatticeDisplayConfig(; kwargs...))
    Compose.set_default_graphic_size(dx, dy)
    return img
end

function _edges(atoms, blockade_radius)
    n = length(atoms)
    edges = Tuple{Int,Int}[]
    for i=1:n, j=1:n
        if sum(abs2, atoms[i] .- atoms[j]) <= blockade_radius ^ 2
            push!(edges, (i, j))
        end
    end
    return edges
end

# Returns a 2-tuple of (image::Context, size)
function img_atoms(al::AtomList; colors, blockade_radius, texts, config)
    @assert length(al) == length(colors)
    @assert length(al) == length(texts)
    atoms = padydim(al).atoms
    rescaler = get_rescaler(atoms, 2*config.scale)
    X = rescaler.xmax - rescaler.xmin + 2*rescaler.pad
    Y = rescaler.ymax - rescaler.ymin + 2*rescaler.pad
    img = _viz_atoms(rescaler.(atoms), _edges(atoms, blockade_radius), colors, texts, config, getscale(rescaler)*blockade_radius)
    img_rescale = 12cm/max(X, Y)
    return Compose.compose(context(0, 0, 1.0, X/Y), img), (X*img_rescale, Y*img_rescale)
end

Base.@kwdef struct LatticeDisplayConfig
    scale::Real = 1.5
    # bond
    bond_color::String = "black"
    # blockade
    blockade_style::String = "none"
end

function _viz_atoms(locs, edges, colors, texts, config, blockade_radius)
    node_styles = [default_node_style(config.scale; color=color) for color in colors]
    edge_style = default_bond_style(config.scale; color=config.bond_color)
    blockade_radius_style = default_blockade_style(config.scale, (config.blockade_style=="half" ? blockade_radius/2 : blockade_radius))
    text_style = default_text_style(config.scale)
    Viznet.canvas() do
        for (i, node) in enumerate(locs)
            node_styles[i] >> node
            if config.blockade_style != "none"
                blockade_radius_style >> node
            end
            text_style >> (node, texts[i])
        end
        for (i, j) in edges
            edge_style >> (locs[i], locs[j])
        end
    end
end

"""
    viz_maskedgrid([io::Union{IO,AbstractString}, ]maskedgrid::MaskedGrid;
        format=PNG,
        blockade_radius = blockade_radius,
        colors=fill("white", count(maskedgrid.mask)))

Draw a `maskedgrid` with colors specified by `colors` and texts specified by `texts`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.
"""
function viz_maskedgrid(io, maskedgrid::MaskedGrid;
        format=PNG,
        colors=fill("white", count(maskedgrid.mask)),
        texts = ["$i" for i=1:sum(maskedgrid.mask)],
        blockade_radius = 0,
        kwargs...
        )
    img, (dx, dy) = img_maskedgrid(maskedgrid; colors=colors, texts=texts, blockade_radius=blockade_radius, config=LatticeDisplayConfig(; kwargs...))
    Compose.draw(format(io, dx, dy), img)
    return
end

# Returns a 2-tuple of (image::Context, size)
function img_maskedgrid(maskedgrid::MaskedGrid; colors, texts, config, blockade_radius)
    atoms = collect_atoms(maskedgrid)
    rescaler = get_rescaler(atoms, 2*config.scale)
    X = rescaler.xmax - rescaler.xmin + rescaler.pad
    Y = rescaler.ymax - rescaler.ymin + rescaler.pad
    line_style_grid = default_line_style_grid(config.scale)
    img1 = _viz_atoms(rescaler.(atoms), _edges(atoms, blockade_radius), colors, texts, config, getscale(rescaler)*blockade_radius)
    img2 = _viz_grid(rescaler.(maskedgrid.xs; dims=1), rescaler.(maskedgrid.ys; dims=2), line_style_grid, Y/X)
    img = compose(context(0, 0, 1.0, X/Y), (context(), img1), (context(), img2))
    img_rescale = 12cm/max(X, Y)
    return img, (X*img_rescale, Y*img_rescale)
end
function _viz_grid(xs, ys, line_style, ymax)
    Viznet.canvas() do
        for i=1:length(xs)
            line_style >> ((xs[i], 0.0), (xs[i], ymax))
        end
        for i=1:length(ys)
            line_style >> ((0.0, ys[i]), (1.0, ys[i]))
        end
    end
end

for (mime, format) in [MIME"image/png"=>PNG, MIME"text/html"=>SVG]
    @eval begin
        function Base.show(io::IO, ::$mime, lt::AbstractLattice{D}) where D
            show(io, $mime(), generate_sites(lt, ntuple(i->5, D)...))
        end
    
        function Base.show(io::IO, ::$mime, maskedgrid::MaskedGrid)
            viz_maskedgrid(io, maskedgrid; format=$format)
        end
    
        function Base.show(io::IO, ::$mime, list::AtomList)
            viz_atoms(io, list; format=$format)
        end
    end
end