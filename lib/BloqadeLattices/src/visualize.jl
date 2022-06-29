const DEFAULT_LINE_COLOR = Ref("#000000")
const DEFAULT_TEXT_COLOR = Ref("#000000")
const DEFAULT_NODE_COLOR = Ref("#FFFFFF")

const CONFIGHELP = """
### Extra Keyword Arguments
* `scale::Float64 = 1.0` is the overall scaling
* `pad::Float64 = 1.0` is the padding space

##### axes
* `axes_text_color::String = DEFAULT_TEXT_COLOR[]`
* `axes_text_fontsize::Float64 = 11.0`
* `axes_num_of_xticks = 5`
* `axes_num_of_yticks = 5`
* `axes_x_offset::Float64 = 0.1`
* `axes_y_offset::Float64 = 0.06`
* `axes_unit::String = "μm"`

##### node
* `node_text_fontsize::Float64 = 5.0`
* `node_text_color::String = DEFAULT_TEXT_COLOR[]`
* `node_stroke_color = DEFAULT_LINE_COLOR[]`
* `node_stroke_linewidth = 0.03`
* `node_fill_color = DEFAULT_NODE_COLOR[]`

##### bond
* `bond_color::String = DEFAULT_LINE_COLOR[]`
* `bond_linewidth::Float64 = 0.03`
# `blockade`
* `blockade_radius::Float64=0`
* `blockade_style::String = "none"`
* `blockade_stroke_color::String = DEFAULT_LINE_COLOR[]`
* `blockade_fill_color::String = "transparent"`
* `blockade_fill_opacity::Float64 = 0.5`
* `blockade_stroke_linewidth = 0.03

##### arrow
* `arrow_linewidth`
* `arrow_color`
* `arrow_head_length`
"""


function config_plotting(sites, pad)
    axes_x_offset = 0.5 * pad
    axes_y_offset = 0.4 * pad
    xmin = minimum(first, sites)
    xmax = maximum(first, sites)
    ymin = minimum(last, sites)
    ymax = maximum(last, sites)
    axes_num_of_yticks = ceil(Int, min((ymax - ymin + 1e-5), 5))
    axes_num_of_xticks = ceil(Int, min((xmax - xmin + 1e-5), 5))
    return (
        axes_x_offset = axes_x_offset,
        axes_y_offset = axes_y_offset,
        axes_num_of_xticks = axes_num_of_xticks,
        axes_num_of_yticks = axes_num_of_yticks,
    )
end

"""
    img_atoms(atoms::AtomList;
        colors=[DEFAULT_LINE_COLOR[], ...],
        texts=["1", "2", ...],
        vectors = [],
        format=:svg,
        filename=nothing,
        kwargs...
        )

Plots `atoms` with colors specified by `colors` and texts specified by `texts`.
Extra vectors can be specified by the `vectors` keyword argument, which is a vector of (start_loc, end_loc) pair.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.
If you want to write this image to the disk without displaying it in a frontend, please try

```julia
julia> img_atoms(generate_sites(SquareLattice(), 5, 5); filename="test.png")
```

The `format` keyword argument can be `:svg`, `:pdf` or `:png`.
Atoms within `blockade_radius` will be connected by edges.

$CONFIGHELP
"""
function img_atoms(
    atoms::AtomList{2};
    colors = nothing,
    texts = nothing,
    vectors = [],
    format = :svg,
    filename = nothing,
    pad=1.0,
    scale=1.0,
    kwargs...,
)
    if length(atoms) == 0
        LuxorGraphPlot._draw(f, 100, 100; format, filename)
    else
        canvas = LuxorGraphPlot.config_canvas(atoms, pad)
        unit = 60 * scale
        Dx, Dy = (canvas.xspan+2*pad)*unit, (canvas.yspan+2*pad)*unit
        config = LatticeDisplayConfig(; pad, config_plotting(atoms, pad)..., kwargs...)
        LuxorGraphPlot._draw(Dx, Dy; format, filename) do
            _viz_atoms(map(loc->(loc[1]+canvas.offsetx, loc[2]+canvas.offsety), atoms), _edges(atoms, config.blockade_radius),
                colors, vectors, texts, config)
            _viz_axes(; canvas..., axes_num_of_xticks=config.axes_num_of_xticks, axes_num_of_yticks=config.axes_num_of_yticks,
                axes_text_fontsize=config.axes_text_fontsize, axes_text_color=config.axes_text_color, axes_unit=config.axes_unit)
        end
    end
end
img_atoms(atoms::AtomList{1}; kwargs...) = img_atoms(padydim(atoms); kwargs...)

function _edges(atoms, blockade_radius)
    n = length(atoms)
    edges = Tuple{Int,Int}[]
    for i in 1:n, j in i+1:n
        if sum(abs2, atoms[i] .- atoms[j]) <= blockade_radius^2
            push!(edges, (i, j))
        end
    end
    return edges
end

_LinRange(x, y, n) = n > 1 ? LinRange(x, y, n) : (x + y) / 2
function _viz_axes(; offsetx, xspan, offsety, yspan, axes_num_of_xticks, axes_num_of_yticks, axes_text_fontsize, axes_text_color, axes_unit)
    xs = _LinRange(offsetx, offsetx+xspan, axes_num_of_xticks)
    ys = _LinRange(offsety, offsety+yspan, axes_num_of_yticks)
    xys = [xs..., ys...]
    locs = [[(x, offsety) for x in xs]..., [(offsety, y) for y in ys]...]
    for (x, loc) in zip(xys, locs)
        LuxorGraphPlot.draw_text(LuxorGraphPlot.Point(loc), "$(round(x; digits=2))$(axes_unit)"; color=axes_text_color, fontsize=axes_text_fontsize)
    end
end

Base.@kwdef struct LatticeDisplayConfig
    # line, node and text
    scale::Float64 = 1.0
    pad::Float64 = 1.0

    # axes
    axes_text_color::String = DEFAULT_LINE_COLOR[]  # NOTE: follow the line color!
    axes_text_fontsize::Float64 = 11.0
    axes_num_of_xticks::Int = 5
    axes_num_of_yticks::Int = 5
    axes_x_offset::Float64 = 0.1
    axes_y_offset::Float64 = 0.06
    axes_unit::String = "μm"

    # node
    node_size::Float64 = 0.15
    node_text_fontsize::Float64 = 5.0
    node_text_color::String = DEFAULT_TEXT_COLOR[]
    node_stroke_color::String = DEFAULT_LINE_COLOR[]
    node_stroke_linewidth::Float64 = 0.03
    node_fill_color::String = DEFAULT_NODE_COLOR[]

    # bond
    bond_color::String = DEFAULT_LINE_COLOR[]
    bond_linewidth::Float64 = 0.03

    # blockade
    blockade_style::String = "none"
    blockade_stroke_color::String = DEFAULT_LINE_COLOR[]
    blockade_fill_color::String = "transparent"
    blockade_fill_opacity::Float64 = 0.5
    blockade_stroke_linewidth::Float64 = 0.03
    blockade_radius::Float64 = 0.0

    # arrow
    arrow_linewidth::Float64=1.0
    arrow_head_length::Float64=6.0
    arrow_color::String="black"
end

function _viz_atoms(locs, edges, colors, vectors, texts, config)
    # show the arrows
    LuxorGraphPlot.Luxor.sethue(config.arrow_color)
    LuxorGraphPlot.Luxor.setline(config.arrow_linewidth)
    for v in vectors
        Luxor.arrow(Point(v[1]), Point(v[2]), arrowheadlength=config.arrow_head_length)
    end

    # show the graph
    graphconfig = LuxorGraphPlot.GraphDisplayConfig(; vertex_stroke_color=config.node_stroke_color, vertex_line_width=config.node_stroke_linewidth,
        vertex_size=config.node_size,
        vertex_text_color=config.node_text_color, fontsize=config.node_text_fontsize,
        edge_color=config.bond_color, edge_line_width=config.bond_linewidth)
    LuxorGraphPlot._show_graph(locs, edges, colors, nothing, nothing, nothing, nothing, nothing, texts, graphconfig)

    # show the blockade
    if config.blockade_style != "none"
        blockadeconfig = LuxorGraphPlot.GraphDisplayConfig(; vertex_stroke_color=config.blockade_stroke_color, vertex_line_width=config.blockade_stroke_linewidth,
        vertex_color=blockade_fill_color, vertex_size=radi, line_style="dashed", fill_opacity=config.blockade_fill_opacity)
        LuxorGraphPlot._show_graph(locs, Tuple{Int,Int}[], nothing, nothing, nothing, nothing, nothing, nothing, texts, blockadeconfig)
    end
end

struct ByDensity
    values::Vector{Float64}
    colormap::String
    vmin::Float64
    vmax::Float64
end

"""
    ByDensity(values; colormap="Grays", vmin=minimum(values), vmax=maximum(values))

For specifying the colors for density plots, where `values` are densities.

# Keyword arguments
* `colormap` is a string for specifying the color map, check the documentation of [`Colors`] package for the detailed description.
* `vmin` and `vmax` are the color range.
"""
function ByDensity(values; colormap = "Grays", vmin = minimum(values), vmax = maximum(values))
    @assert vmax >= vmin
    return ByDensity(values, colormap, vmin, vmax)
end

function resolve_colors(::Nothing, locs, config)
    return fill(config.node_fill_color, length(locs))
end
function resolve_colors(colors::String, locs, config)
    return fill(colors, length(locs))
end
function resolve_colors(colors, locs, config)
    @assert length(locs) == length(colors)
    return collect(String, colors)
end
function resolve_colors(colors::ByDensity, locs, config)
    @assert length(locs) == length(colors.values)
    N = 100
    cmap = Compose.colormap(colors.colormap, N)
    return map(colors.values) do v
        scale = max(colors.vmax - colors.vmin, 1e-12)  # avoid zero devision
        index = max(min(ceil(Int, (v - colors.vmin) / scale * N), N), 1)
        return cmap[index]
    end
end

"""
    img_maskedgrid(maskedgrid::MaskedGrid;
        format=:svg,
        filename=nothing,
        colors=[DEFAULT_LINE_COLOR[], ...],
        texts=["1", "2", ...],
        vectors=[],
        blockade_radius = 0,
        kwargs...
        )

Draw a `maskedgrid` with colors specified by `colors` and texts specified by `texts`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.

See also the docstring of [`img_atoms`](@ref) for explanations of other keyword arguments.

$CONFIGHELP
"""
function img_maskedgrid(
    maskedgrid::MaskedGrid;
    format = :svg,
    filename = nothing,
    colors = nothing,
    texts = nothing,
    vectors = [],
    blockade_radius = 0,
    kwargs...,
)
    atoms = padydim(collect_atoms(maskedgrid))
    isempty(atoms) && return
    img, (dx, dy) = viz_maskedgrid(
        maskedgrid;
        colors,
        vectors,
        texts,
        blockade_radius,
        config = LatticeDisplayConfig(; config_plotting(atoms, pad)..., kwargs...),
    )
    if length(locations) == 0
        _draw(f, 100, 100; format, filename)
    else
        config = LuxorGraphPlot.autoconfig(locations; pad, kwargs...)
        Dx, Dy = (config.xspan+2*config.pad)*config.unit, (config.yspan+2*config.pad)*config.unit
        _draw(Dx, Dy; format, filename) do
            _show_graph(map(loc->(loc[1]+config.offsetx, loc[2]+config.offsety), locations), edges,
            vertex_colors, vertex_stroke_colors, vertex_text_colors, vertex_sizes, vertex_shapes, edge_colors, texts, config)
            f()
        end
    end
end

# Returns a 2-tuple of (image::Context, size)
function viz_maskedgrid(maskedgrid::MaskedGrid; colors, vectors, texts, config, blockade_radius)
    atoms = padydim(collect_atoms(maskedgrid))
    rescaler = get_rescaler(atoms, config.pad)
    rescale = getscale(rescaler) * config.image_size * config.scale * 1.6
    line_style_grid = Viznet.bondstyle(:default, Compose.stroke("#AAAAAA"), linewidth(0.3mm * rescale); dashed = true)
    img1 = _viz_atoms(
        rescaler.(atoms),
        _edges(atoms, blockade_radius),
        colors,
        vectors,
        texts,
        config,
        blockade_radius,
        getscale(rescaler),
    )
    ymax = (rescaler.ymax - rescaler.ymin + 2 * rescaler.pad) / (rescaler.xmax - rescaler.xmin + 2 * rescaler.pad)
    img2 = _viz_grid(rescaler.(maskedgrid.xs; dims = 1), rescaler.(maskedgrid.ys; dims = 2), line_style_grid, ymax)
    img_axes = _viz_axes(rescaler, config)
    return fit_image(rescaler, config.image_size, img1, img2, img_axes)
end
function _viz_grid(xs, ys, line_style, ymax)
    Viznet.canvas() do
        for i in 1:length(xs)
            line_style >> ((xs[i], 0.0), (xs[i], ymax))
        end
        for i in 1:length(ys)
            line_style >> ((0.0, ys[i]), (1.0, ys[i]))
        end
    end
end

for (mime, format) in [MIME"image/png" => :png, MIME"text/html" => :svg]
    @eval begin
        function Base.show(io::IO, m::$mime, maskedgrid::MaskedGrid)
            Base.show(io, m, img_maskedgrid(maskedgrid; format = $(QuoteNode(format))))
            return nothing
        end

        function Base.show(io::IO, m::$mime, list::AtomList)
            Base.show(io, m, img_atoms(list; format = $(QuoteNode(format))))
            return nothing
        end
    end
end

function darktheme!()
    DEFAULT_LINE_COLOR[] = "#FFFFFF"
    DEFAULT_TEXT_COLOR[] = "#FFFFFF"
    DEFAULT_NODE_COLOR[] = "#000000"
end

function lighttheme!()
    DEFAULT_LINE_COLOR[] = "#000000"
    DEFAULT_TEXT_COLOR[] = "#000000"
    DEFAULT_NODE_COLOR[] = "#FFFFFF"
end