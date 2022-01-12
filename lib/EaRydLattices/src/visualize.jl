struct Rescaler{T}
    xmin::T
    xmax::T
    ymin::T
    ymax::T
end

getscale(r::Rescaler) = 1/(r.xmax-r.xmin+1)

function (r::Rescaler{T})(x; dims=(1,2)) where T
    xmin, ymin, xmax, ymax = r.xmin, r.ymin, r.xmax, r.ymax
    scale = getscale(r)
    if dims == (1,2)
        return (x[1]-xmin+0.5, ymax+0.5-x[2]) .* scale
    elseif dims == 1
        return (x - xmin + 0.5) * scale
    elseif dims == 2
        return (ymax + 0.5 - x) * scale
    else
        throw(ArgumentError("dims should be (1,2), 1 or 2."))
    end
end

function get_rescaler(atoms::Vector{<:Tuple})
    xmin = minimum(x->x[1], atoms)
    ymin = minimum(x->x[2], atoms)
    xmax = maximum(x->x[1], atoms)
    ymax = maximum(x->x[2], atoms)
    return Rescaler(xmin, xmax, ymin, ymax)
end

default_node_style(scale) = compose(context(), Viznet.nodestyle(:default, r=0.15cm*scale), Compose.stroke("black"), fill("white"), linewidth(0.3mm*scale))
default_text_style(scale) = Viznet.textstyle(:default, fontsize(4pt*scale))
default_line_style_grid(scale) = Viznet.bondstyle(:default, Compose.stroke("#AAAAAA"), linewidth(0.3mm*scale); dashed=true)

"""
    viz_atoms(io::Union{IO,AbstractString}, atoms::Vector{<:Tuple}; scale=1.0)

Draw `atoms` with scaling factor `scale`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.
If you want to write this image to the disk without using a frontend, please check [`img_atoms`](@ref).
"""
function viz_atoms(io, atoms::Vector{<:Tuple}; scale=1.0, format=PNG)
    img, (dx, dy) = img_atoms(atoms; scale=scale)
    Compose.draw(format(io, dx, dy), img)
end

# Returns a 2-tuple of (image::Context, size)
function img_atoms(atoms::Vector{<:Tuple}; scale=1.0)
    rescaler = get_rescaler(atoms)
    xspan = rescaler.xmax - rescaler.xmin
    yspan = rescaler.ymax - rescaler.ymin
    X = (xspan+1)
    Y = (yspan+1)
    node_style = default_node_style(scale)
    text_style = default_text_style(scale)
    img = _viz_atoms(rescaler.(atoms), node_style, text_style)
    return Compose.compose(context(0, 0, 1.0, X/Y), img), (X*scale*cm, Y*scale*cm)
end
function _viz_atoms(locs, node_style, text_style)
    Viznet.canvas() do
        for (i, node) in enumerate(locs)
            node_style >> node
            text_style >> (node, "$i")
        end
    end
end

"""
    viz_maskedgrid(io::Union{IO,AbstractString}, maskedgrid::MaskedGrid; scale=1.0)

Draw a `maskedgrid` with scaling factor `scale`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.
If you want to write this image to the disk without using a frontend, please check [`img_atoms`](@ref).
"""
function viz_maskedgrid(io, maskedgrid::MaskedGrid; scale=1.0, format=PNG)
    img, (dx, dy) = img_maskedgrid(maskedgrid; scale=scale)
    Compose.draw(format(io, dx, dy), img)
end

# Returns a 2-tuple of (image::Context, size)
function img_maskedgrid(maskedgrid::MaskedGrid; scale=1.0)
    atoms = locations(maskedgrid)
    rescaler = get_rescaler(atoms)
    xspan = rescaler.xmax - rescaler.xmin
    yspan = rescaler.ymax - rescaler.ymin
    X = (xspan+1)
    Y = (yspan+1)
    node_style = default_node_style(scale)
    text_style = default_text_style(scale)
    line_style_grid = default_line_style_grid(scale)
    img1 = _viz_atoms(rescaler.(atoms), node_style, text_style)
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
        function Base.show(io::IO, ::$mime, lt::AbstractLattice{2})
            viz_atoms(io, generate_sites(lt, 5, 5); scale=2.0, format=$format)
        end
    
        function Base.show(io::IO, ::$mime, lt::AbstractLattice{1})
            viz_atoms(io, padydim.(generate_sites(lt, 5)); scale=2.0, format=$format)
        end
        function Base.show(io::IO, ::$mime, maskedgrid::MaskedGrid)
            viz_maskedgrid(io, maskedgrid; scale=2.0, format=$format)
        end
    
        function Base.show(io::IO, ::$mime, list::AtomList{Tuple{<:Real}})
            viz_atoms(io, padydim.(list.atoms); scale=2.0, format=$format)
        end
    
        function Base.show(io::IO, ::$mime, list::AtomList)
            viz_atoms(io, list.atoms; scale=2.0, format=$format)
        end
    end
end
