const DEFAULT_LINE_COLOR = Ref("#000000")
const DEFAULT_TEXT_COLOR = Ref("#000000")
const DEFAULT_NODE_COLOR = Ref("transparent")

const CONFIGHELP = """
### Extra Keyword Arguments
* `scale::Float64 = 1.0` is the overall scaling
* `xpad::Float64 = 2.5` is the padding space in x axis
* `ypad::Float64 = 1.5` is the padding space in y axis

##### axes
* `axes_text_color = DEFAULT_TEXT_COLOR[]`
* `axes_text_fontsize::Float64 = 18.0`
* `axes_num_of_xticks = 5`
* `axes_num_of_yticks = 5`
* `axes_x_offset::Float64 = 0.1`
* `axes_y_offset::Float64 = 0.06`
* `axes_unit::String = "μm"`

##### node
* `node_text_fontsize::Float64 = 14.0`
* `node_text_color = DEFAULT_TEXT_COLOR[]`
* `node_stroke_color = DEFAULT_LINE_COLOR[]`
* `node_stroke_linewidth = 1`
* `node_fill_color = DEFAULT_NODE_COLOR[]`

##### bond
* `bond_color = DEFAULT_LINE_COLOR[]`
* `bond_linewidth::Float64 = 1.0`
# `blockade`
* `blockade_radius::Float64=0`
* `blockade_style::String = "none"`
* `blockade_stroke_color = DEFAULT_LINE_COLOR[]`
* `blockade_fill_color = "transparent"`
* `blockade_fill_opacity::Float64 = 0.5`
* `blockade_stroke_linewidth = 1.0`   # in pt

##### arrow
* `arrow_linewidth`
* `arrow_color`
* `arrow_head_length`

##### grid
* grid_stroke_color="#AAAAAA"
* grid_stroke_width::Float64=1
* grid_stroke_style::String="dashed"
"""


function config_plotting(sites, xpad, ypad)
    axes_x_offset = 0.5 * xpad
    axes_y_offset = 0.4 * ypad
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
        xpad=1.5,
        ypad=1.0,
        scale=1.0,
        kwargs...,
    )
    length(atoms) == 0 && return LuxorGraphPlot._draw(()->nothing, 100, 100; format, filename)

    canvas = LuxorGraphPlot.config_canvas(atoms, xpad, ypad)
    xmin = minimum(first, atoms)
    ymin = minimum(last, atoms)
    unit = ceil(Int, 60 * scale)
    Dx, Dy = (canvas.xspan+2*xpad)*unit, (canvas.yspan+2*ypad)*unit
    config = LatticeDisplayConfig(; scale, xpad, ypad, config_plotting(atoms, xpad, ypad)..., kwargs...)
    transform(loc) = (loc[1]+canvas.offsetx, loc[2]+canvas.offsety)
    LuxorGraphPlot._draw(Dx, Dy; format, filename) do
        _viz_atoms(transform.(atoms), _edges(atoms, config.blockade_radius),
            colors, map(ab->transform.(ab), vectors), texts, config)
        _viz_axes(; xmin, ymin, axes_x_offset=config.axes_x_offset, axes_y_offset=config.axes_y_offset, canvas...,
            axes_num_of_xticks=config.axes_num_of_xticks, axes_num_of_yticks=config.axes_num_of_yticks,
            axes_text_fontsize=config.axes_text_fontsize, axes_text_color=config.axes_text_color, axes_unit=config.axes_unit, scale=scale)
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
function _viz_axes(; offsetx, xmin, xspan, offsety, ymin, yspan,
        axes_num_of_xticks, axes_num_of_yticks, axes_text_fontsize, axes_text_color, axes_unit, scale,
        axes_x_offset, axes_y_offset
        )
    unit = 60 * scale
    # the true coordinates
    xs = _LinRange(xmin, xmin+xspan, axes_num_of_xticks)
    ys = _LinRange(ymin, ymin+yspan, axes_num_of_yticks)
    xys = [xs..., ys...]
    # the display coordinates of axes labels
    locs = [[Point(offsetx+x, ymin+yspan+offsety+axes_y_offset)*unit for x in xs]..., [Point(offsetx+xmin-axes_x_offset, offsety+y)*unit for y in ys]...]
    for (x, loc) in zip(xys, locs)
        LuxorGraphPlot.draw_text(loc, "$(round(x; digits=2))$(axes_unit)"; color=axes_text_color, fontsize=axes_text_fontsize)
    end
end

Base.@kwdef struct LatticeDisplayConfig
    # line, node and text
    scale::Float64 = 1.0
    xpad::Float64 = 2.5
    ypad::Float64 = 1.5

    # axes
    axes_text_color = DEFAULT_LINE_COLOR[]  # NOTE: follow the line color!
    axes_text_fontsize::Float64 = 18.0
    axes_num_of_xticks::Int = 5
    axes_num_of_yticks::Int = 5
    axes_x_offset::Float64 = 0.1
    axes_y_offset::Float64 = 0.06
    axes_unit::String = "μm"

    # node
    node_size::Float64 = 0.25
    node_text_fontsize::Float64 = 14.0
    node_text_color = DEFAULT_TEXT_COLOR[]
    node_stroke_color = DEFAULT_LINE_COLOR[]
    node_stroke_linewidth::Float64 = 1.0   # in pt
    node_fill_color = DEFAULT_NODE_COLOR[]

    # bond
    bond_color = DEFAULT_LINE_COLOR[]
    bond_linewidth::Float64 = 1.0  # in pt

    # blockade
    blockade_style::String = "none"
    blockade_stroke_color = DEFAULT_LINE_COLOR[]
    blockade_fill_color = "transparent"
    blockade_fill_opacity::Float64 = 0.5
    blockade_stroke_linewidth::Float64 = 1.0
    blockade_radius::Float64 = 0.0

    # arrow
    arrow_linewidth::Float64=1.0   # in pt
    arrow_head_length::Float64=6.0
    arrow_color="black"

    # grid
    grid_stroke_color="#AAAAAA"
    grid_stroke_width::Float64=1
    grid_stroke_style::String="dashed"
end

function _viz_atoms(locs, edges, colors, vectors, texts, config)
    # show the blockade
    if config.blockade_style != "none"
        blockadeconfig = LuxorGraphPlot.GraphDisplayConfig(; vertex_stroke_color=config.blockade_stroke_color, vertex_line_width=config.blockade_stroke_linewidth,
            vertex_fill_color=Colors.RGBA(LuxorGraphPlot.sethue(config.blockade_fill_color)..., config.blockade_fill_opacity),
            vertex_size=config.blockade_style=="half" ? config.blockade_radius/2 : config.blockade_radius,
            vertex_line_style="dashed",
            #vertex_stroke_opacity=config.blockade_fill_opacity,
            unit=config.scale * 60)
        LuxorGraphPlot._show_graph(locs, Tuple{Int,Int}[], nothing, nothing, nothing, nothing, nothing, nothing, texts, blockadeconfig)
    end

    # show the arrows
    unit = 60 * config.scale
    for v in vectors
        LuxorGraphPlot.draw_edge(Point(v[1])*unit, Point(v[2])*unit; color=config.arrow_color, line_width=config.arrow_linewidth,
            arrow=true, arrowheadlength=config.arrow_head_length, line_style="solid")
    end

    # show the graph
    graphconfig = LuxorGraphPlot.GraphDisplayConfig(; vertex_stroke_color=config.node_stroke_color, vertex_line_width=config.node_stroke_linewidth,
        vertex_size=config.node_size, vertex_fill_color=config.node_fill_color,
        vertex_text_color=config.node_text_color, fontsize=config.node_text_fontsize,
        edge_color=config.bond_color, edge_line_width=config.bond_linewidth, unit=config.scale * 60)
    LuxorGraphPlot._show_graph(locs, edges, colors, nothing, nothing, nothing, nothing, nothing, texts, graphconfig)
end

struct ByDensity
    values::Vector{Float64}
    colormap
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
    return ByDensity(values, Colors.colormap(colormap), vmin, vmax)
end

function LuxorGraphPlot._get(colors::ByDensity, i::Int, default)
    N = 100
    scale = max(colors.vmax - colors.vmin, 1e-12)  # avoid zero devision
    index = max(min(ceil(Int, (colors.values[i] - colors.vmin) / scale * N), N), 1)
    return colors.colormap[index]
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
        scale = 1.0,
        xpad = 2.5,
        ypad = 1.5,
        kwargs...,
    )
    atoms = padydim(collect_atoms(maskedgrid))
    canvas = LuxorGraphPlot.config_canvas(atoms, xpad, ypad)
    xmin = minimum(first, atoms)
    ymin = minimum(last, atoms)
    unit = ceil(Int, 60 * scale)

    Dx, Dy = (canvas.xspan+2*xpad)*unit, (canvas.yspan+2*ypad)*unit
    config = LatticeDisplayConfig(; scale, xpad, ypad, config_plotting(atoms, xpad, ypad)..., kwargs...)
    LuxorGraphPlot._draw(Dx, Dy; format, filename) do
        # show the grid
        _viz_grid(maskedgrid.xs, maskedgrid.ys; scale, xpad, ypad, xmin, ymin, canvas...,
            color=config.grid_stroke_color, line_width=config.grid_stroke_width, line_style=config.grid_stroke_style)
        # show atoms
        length(atoms) > 0 && _viz_atoms(map(loc->(loc[1]+canvas.offsetx, loc[2]+canvas.offsety), atoms), _edges(atoms, config.blockade_radius),
                colors, vectors, texts, config)
        # show the axes
        _viz_axes(; xmin, ymin, axes_x_offset=config.axes_x_offset, axes_y_offset=config.axes_y_offset, canvas...,
            axes_num_of_xticks=config.axes_num_of_xticks, axes_num_of_yticks=config.axes_num_of_yticks,
            axes_text_fontsize=config.axes_text_fontsize, axes_text_color=config.axes_text_color, axes_unit=config.axes_unit, scale=scale)
    end
end

function _viz_grid(xs, ys; xpad, ypad, scale, xmin, ymin, offsetx, offsety, xspan, yspan, color, line_width, line_style)
    _xmax = xmin + xspan + offsetx + 0.3*xpad
    _xmin = xmin + offsetx - 0.3*xpad
    _ymax = ymin + yspan + offsety + 0.3*ypad
    _ymin = ymin + offsety - 0.3*ypad
    unit = scale * 60
    for x in xs
        LuxorGraphPlot.draw_edge(Point(x+offsetx, _ymin)*unit, Point(x+offsetx, _ymax)*unit; color, line_width, line_style)
    end
    for y in ys
        LuxorGraphPlot.draw_edge(Point(_xmin, y+offsety)*unit, Point(_xmax, y+offsety)*unit; color, line_width, line_style)
    end
end

for (mime, format) in [MIME"image/png" => :png, MIME"image/svg+xml" => :svg]
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
    #DEFAULT_NODE_COLOR[] = "#000000"
end

function lighttheme!()
    DEFAULT_LINE_COLOR[] = "#000000"
    DEFAULT_TEXT_COLOR[] = "#000000"
    #DEFAULT_NODE_COLOR[] = "#FFFFFF"
end