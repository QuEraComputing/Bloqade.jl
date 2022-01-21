struct Rescaler{T}
    xmin::T
    xmax::T
    ymin::T
    ymax::T
    pad::T
end

getscale(r::Rescaler) = min(1/(r.xmax-r.xmin+2*r.pad), 1/(r.ymax-r.ymin+2*r.pad))

function (r::Rescaler{T})(x; dims=(1,2)) where T
    xmin, ymin, xmax, ymax, pad = r.xmin, r.ymin, r.xmax, r.ymax, r.pad
    scale = getscale(r)
    if dims == (1,2)
        return (x[1]-xmin+pad, ymax+pad-x[2]) .* scale
    elseif dims == 1
        return (x - xmin + pad) * scale
    elseif dims == 2
        return (ymax + pad - x) * scale
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

default_node_style(scale, stroke_color, fill_color) = compose(context(), Viznet.nodestyle(:default, r=0.15cm*scale), Compose.stroke(stroke_color), fill(fill_color), linewidth(0.3mm*scale))
default_text_style(scale, color) = Viznet.textstyle(:default, fontsize(4pt*scale), fill(color))
default_bond_style(scale, color) = Viznet.bondstyle(:default, Compose.stroke(color), linewidth(0.3mm*scale))
default_line_style_grid(scale) = Viznet.bondstyle(:default, Compose.stroke("#AAAAAA"), linewidth(0.3mm*scale); dashed=true)

"""
    img_atoms(atoms::AtomList;
        colors=["black", "black", ...], blockade_radius=0,
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
function img_atoms(atoms::AtomList;
        colors=nothing,
        blockade_radius=0,
        texts = nothing,
        format=SVG, io=nothing,
        kwargs...)
    img, (dx, dy) = viz_atoms(atoms; colors=colors, blockade_radius=blockade_radius, texts=texts, config=LatticeDisplayConfig(; kwargs...))
    if io === nothing
        Compose.set_default_graphic_size(dx, dy)
        return img
    else
        return format(io, dx, dy)(img)
    end
end

function _edges(atoms, blockade_radius)
    n = length(atoms)
    edges = Tuple{Int,Int}[]
    for i=1:n, j=i+1:n
        if sum(abs2, atoms[i] .- atoms[j]) <= blockade_radius ^ 2
            push!(edges, (i, j))
        end
    end
    return edges
end

function fit_image(rescaler::Rescaler, image_size, imgs...)
    X = rescaler.xmax - rescaler.xmin + 2*rescaler.pad
    Y = rescaler.ymax - rescaler.ymin + 2*rescaler.pad
    img_rescale = image_size/max(X, Y)*cm
    if Y < X
        return Compose.compose(context(0, 0, 1.0, X/Y), imgs...), (X*img_rescale, Y*img_rescale)
    else
        return Compose.compose(context(0, 0, Y/X, 1.0), imgs...), (X*img_rescale, Y*img_rescale)
    end
end

# Returns a 2-tuple of (image::Context, size)
function viz_atoms(al::AtomList; colors, blockade_radius, texts, config)
    atoms = padydim(al).atoms
    rescaler = get_rescaler(atoms, config.pad)
    img = _viz_atoms(rescaler.(atoms), _edges(atoms, blockade_radius), colors, texts, config, blockade_radius, getscale(rescaler))
    img_axes = _viz_axes(rescaler, config)
    return fit_image(rescaler, config.image_size, img, img_axes)
end

function _viz_axes(rescaler, config)
    xs = LinRange(rescaler.xmin, rescaler.xmax, config.axes_num_of_xticks)
    ys = LinRange(rescaler.ymin, rescaler.ymax, config.axes_num_of_yticks)
    xlocs = rescaler.([(x, -config.axes_y_offset) for x in xs])
    ylocs = rescaler.([(-config.axes_x_offset, y) for y in ys])
    return _axes!([xs..., ys...], [xlocs..., ylocs...], config, getscale(rescaler))
end

Base.@kwdef struct LatticeDisplayConfig
    # line, node and text
    scale::Float64 = 1.0
    pad::Float64 = 1.5

    # axes
    axes_text_color::String = "black"
    axes_text_fontsize::Float64 = 6.0
    axes_num_of_xticks = 5
    axes_num_of_yticks = 5
    axes_x_offset::Float64 = 0.75
    axes_y_offset::Float64 = 0.5
    axes_unit::String = "μm"

    # node
    node_text_color::String = "black"
    node_stroke_color = "black"
    node_fill_color = "white"
    # bond
    bond_color::String = "black"
    # blockade
    blockade_style::String = "none"
    blockade_stroke_color::String = "black"
    blockade_fill_color::String = "transparent"
    blockade_fill_opacity::Float64 = 0.5
    # image size in cm
    image_size::Float64 = 12
end

function _viz_atoms(locs, edges, colors, texts, config, blockade_radius, rescale)
    radi = (config.blockade_style=="half" ? blockade_radius/2 : blockade_radius)*rescale
    rescale = rescale * config.image_size * config.scale * 1.6
    if colors !== nothing
        @assert length(locs) == length(colors)
        node_styles = [default_node_style(rescale, config.node_stroke_color, color) for color in colors]
    else
        node_styles = fill(default_node_style(rescale, config.node_stroke_color, config.node_fill_color), length(locs))
    end
    if texts !== nothing
        @assert length(locs) == length(texts)
    end
    edge_style = default_bond_style(rescale, config.bond_color)
    blockade_radius_style = Viznet.nodestyle(:circle,
        Compose.stroke(config.blockade_stroke_color),
        Compose.strokedash([0.5mm*rescale, 0.5mm*rescale]),
        Compose.linewidth(0.2mm*rescale),
        Compose.fill(@show config.blockade_fill_color), 
        Compose.fillopacity(config.blockade_fill_opacity);
        r=radi)
    text_style = default_text_style(rescale, config.node_text_color)
    Viznet.canvas() do
        for (i, node) in enumerate(locs)
            node_styles[i] >> node
            if config.blockade_style != "none"
                blockade_radius_style >> node
            end
            if config.node_text_color !== "transparent"
                text_style >> (node, texts === nothing ? "$i" : texts[i])
            end
        end
        for (i, j) in edges
            edge_style >> (locs[i], locs[j])
        end
    end
end

function _axes!(xs, locs, config, rescale)
    rescale = rescale * config.image_size * config.scale * 1.6
    text_style = Viznet.textstyle(:default, fontsize((config.axes_text_fontsize)*pt), fill(config.axes_text_color))
    Viznet.canvas() do
        for (x, loc) in zip(xs, locs)
            text_style >> (loc, "$(round(x; digits=2))$(config.axes_unit)")
        end
    end
end

"""
    img_maskedgrid([io::Union{IO,AbstractString}, ]maskedgrid::MaskedGrid;
        format=PNG,
        blockade_radius = blockade_radius,
        colors=fill("white", count(maskedgrid.mask)))

Draw a `maskedgrid` with colors specified by `colors` and texts specified by `texts`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.
"""
function img_maskedgrid(maskedgrid::MaskedGrid;
        format=SVG, io=nothing,
        colors=nothing,
        texts = nothing,
        blockade_radius = 0,
        kwargs...
        )
    img, (dx, dy) = viz_maskedgrid(maskedgrid; colors=colors, texts=texts, blockade_radius=blockade_radius, config=LatticeDisplayConfig(; kwargs...))
    if io === nothing
        Compose.set_default_graphic_size(dx, dy)
        return img
    else
        return format(io, dx, dy)(img)
    end
end

# Returns a 2-tuple of (image::Context, size)
function viz_maskedgrid(maskedgrid::MaskedGrid; colors, texts, config, blockade_radius)
    atoms = collect_atoms(maskedgrid)
    rescaler = get_rescaler(atoms, config.pad)
    line_style_grid = default_line_style_grid(config.scale)
    img1 = _viz_atoms(rescaler.(atoms), _edges(atoms, blockade_radius), colors, texts, config, blockade_radius, getscale(rescaler))
    ymax = (rescaler.ymax - rescaler.ymin + 2*rescaler.pad)/(rescaler.xmax - rescaler.xmin + 2*rescaler.pad)
    img2 = _viz_grid(rescaler.(maskedgrid.xs; dims=1), rescaler.(maskedgrid.ys; dims=2), line_style_grid, ymax)
    img_axes = _viz_axes(rescaler, config)
    fit_image(rescaler, config.image_size, img1, img2, img_axes)
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
            img_maskedgrid(maskedgrid; format=$format, io=io)
            nothing
        end
    
        function Base.show(io::IO, ::$mime, list::AtomList)
            img_atoms(list; format=$format, io=io)
            nothing
        end
    end
end