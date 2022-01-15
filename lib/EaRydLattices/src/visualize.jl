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
default_line_style_grid(scale) = Viznet.bondstyle(:default, Compose.stroke("#AAAAAA"), linewidth(0.3mm*scale); dashed=true)

"""
    viz_atoms(io::Union{IO,AbstractString}, atoms::AtomList; scale=1.0, colors=fill("white", count(maskedgrid.mask)))

Draw `atoms` with scaling factor `scale`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.
If you want to write this image to the disk without using a frontend, please check [`img_atoms`](@ref).
"""
function viz_atoms(io, atoms::AtomList; scale=1.0, format=PNG, colors=fill("white", length(atoms)))
    img, (dx, dy) = img_atoms(atoms; scale=scale, colors=colors)
    @show dx, dy
    Compose.draw(format(io, dx, dy), img)
    return
end

# Returns a 2-tuple of (image::Context, size)
function img_atoms(al::AtomList; scale, colors)
    atoms = padydim(al).atoms
    rescaler = get_rescaler(atoms, 2*scale)
    X = rescaler.xmax - rescaler.xmin + 2*rescaler.pad
    Y = rescaler.ymax - rescaler.ymin + 2*rescaler.pad
    node_styles = [default_node_style(scale; color=color) for color in colors]
    text_style = default_text_style(scale)
    img = _viz_atoms(rescaler.(atoms), node_styles, text_style)
    img_rescale = 12cm/max(X, Y)
    return Compose.compose(context(0, 0, 1.0, X/Y), img), (X*img_rescale, Y*img_rescale)
end
function _viz_atoms(locs, node_styles, text_style)
    Viznet.canvas() do
        for (i, node) in enumerate(locs)
            node_styles[i] >> node
            text_style >> (node, "$i")
        end
    end
end

"""
    viz_maskedgrid(io::Union{IO,AbstractString}, maskedgrid::MaskedGrid; scale=1.0, colors=fill("white", count(maskedgrid.mask)))

Draw a `maskedgrid` with scaling factor `scale`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.
If you want to write this image to the disk without using a frontend, please check [`img_atoms`](@ref).
"""
function viz_maskedgrid(io, maskedgrid::MaskedGrid; scale=1.0, format=PNG, colors=fill("white", count(maskedgrid.mask)))
    img, (dx, dy) = img_maskedgrid(maskedgrid; scale=scale, colors=colors)
    Compose.draw(format(io, dx, dy), img)
    return
end

# Returns a 2-tuple of (image::Context, size)
function img_maskedgrid(maskedgrid::MaskedGrid; scale, colors)
    atoms = collect_atoms(maskedgrid)
    rescaler = get_rescaler(atoms, 2*scale)
    X = rescaler.xmax - rescaler.xmin + rescaler.pad
    Y = rescaler.ymax - rescaler.ymin + rescaler.pad
    node_styles = [default_node_style(scale; color=color) for color in colors]
    text_style = default_text_style(scale)
    line_style_grid = default_line_style_grid(scale)
    img1 = _viz_atoms(rescaler.(atoms), node_styles, text_style)
    img2 = _viz_grid(rescaler.(maskedgrid.xs; dims=1), rescaler.(maskedgrid.ys; dims=2), line_style_grid, Y/X)
    img = compose(context(0, 0, 1.0, X/Y), (context(), img1), (context(), img2))
    return img, (X*scale*cm, Y*scale*cm)
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
            viz_maskedgrid(io, maskedgrid; scale=scale_heuristic(collect_atoms(maskedgrid)), format=$format)
        end
    
        function Base.show(io::IO, ::$mime, list::AtomList)
            viz_atoms(io, list; scale=scale_heuristic(list), format=$format)
        end
    end
end

# get the atom size in the plot
function scale_heuristic(list::AtomList)
    ds = average_distance(list)
    return ds / sqrt(length(list)) * 0.8
end

# get the average distance between atoms
function average_distance(list::AtomList)
    n = length(list)
    n <= 1 && return 1.0
    ds = 0.0
    @inbounds for i=1:n
        for j=i+1:n
            ds += sqrt(sum(abs2, list[i] .- list[j]))
        end
    end
    return ds/(n*(n-1)/2)
end