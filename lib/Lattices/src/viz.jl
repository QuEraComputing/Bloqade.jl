import Viznet
using Viznet.Compose

export viz_atoms, viz_maskedgrid

struct Rescaler{T}
    xmin::T
    xmax::T
    ymin::T
    ymax::T
end

getscale(r::Rescaler) = 1/(max(r.xmax-r.xmin+1, r.ymax-r.ymin+1))

function (r::Rescaler{T})(x; dims=(1,2)) where T
    xmin, ymin = r.xmin, r.ymin
    scale = getscale(r)
    if dims == (1,2)
        return (x .- (xmin-0.5, ymin-0.5)) .* scale
    elseif dims == 1
        return (x - xmin + 0.5) * scale
    elseif dims == 2
        return (x - ymin + 0.5) * scale
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

default_node_style(scale) = compose(context(), Viznet.nodestyle(:default, r=0.15cm*scale), stroke("black"), fill("white"), linewidth(0.3mm*scale))
default_text_style(scale) = Viznet.textstyle(:default, fontsize(4pt*scale))
default_line_style_grid(scale) = Viznet.bondstyle(:default, stroke("#AAAAAA"), linewidth(0.3mm*scale); dashed=true)
function viz_atoms(atoms::Vector{<:Tuple}; scale=1.0)
    img, dx, dy = img_atoms(atoms; scale=scale)
    Compose.draw(SVG(dx, dy), img)
end
function img_atoms(atoms::Vector{<:Tuple}; scale=1.0)
    rescaler = get_rescaler(atoms)
    xspan = rescaler.xmax - rescaler.xmin
    yspan = rescaler.ymax - rescaler.ymin
    X = (xspan+1)
    Y = (yspan+1)
    node_style = default_node_style(scale)
    text_style = default_text_style(scale)
    img = _viz_atoms(rescaler.(atoms), node_style, text_style)
    return Compose.compose(context(0, 0, 1.0, X/Y), img), X*scale*cm, Y*scale*cm
end
function _viz_atoms(locs, node_style, text_style)
    Viznet.canvas() do
        for (i, node) in enumerate(locs)
            node_style >> node
            text_style >> (node, "$i")
        end
    end
end

function viz_maskedgrid(mg::MaskedGrid; scale=1.0)
    img, dx, dy = img_maskedgrid(mg; scale=scale)
    Compose.draw(SVG(dx, dy), img)
end
function img_maskedgrid(mg::MaskedGrid; scale=1.0)
    atoms = locations(mg)
    rescaler = get_rescaler(atoms)
    xspan = rescaler.xmax - rescaler.xmin
    yspan = rescaler.ymax - rescaler.ymin
    X = (xspan+1)
    Y = (yspan+1)
    node_style = default_node_style(scale)
    text_style = default_text_style(scale)
    line_style_grid = default_line_style_grid(scale)
    img1 = _viz_atoms(rescaler.(atoms), node_style, text_style)
    img2 = _viz_grid(rescaler.(mg.xs; dims=1), rescaler.(mg.ys; dims=2), line_style_grid)
    img = compose(context(0, 0, 1.0, X/Y), (context(), img1), (context(), img2))
    return img, X*scale*cm, Y*scale*cm
end
function _viz_grid(xs, ys, line_style)
    Viznet.canvas() do
        for i=1:length(xs)
            line_style >> ((xs[i], 0.0), (xs[i], 1.0))
        end
        for i=1:length(ys)
            line_style >> ((0.0, ys[i]), (1.0, ys[i]))
        end
    end
end

# for Pluto and vscode
for mime in [:(MIME"text/html"), :(MIME"image/png")]
    @eval function Base.show(io::IO, ::$mime, lt::AbstractLattice{2})
        Base.show(io, $mime(), viz_atoms(generate_sites(lt, 5, 5); scale=2.0))
    end
    @eval function Base.show(io::IO, ::$mime, lt::AbstractLattice{1})
        Base.show(io, $mime(), viz_atoms(padydim.(generate_sites(lt, 5)); scale=2.0))
    end
end
padydim(x::Tuple{T}) where T = (x[1], zero(T))