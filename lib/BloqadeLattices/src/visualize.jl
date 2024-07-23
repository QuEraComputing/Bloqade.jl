const DEFAULT_LINE_COLOR = Ref("#000000")
const DEFAULT_TEXT_COLOR = Ref("#000000")
const DEFAULT_NODE_COLOR = Ref("transparent")
const DEFAULT_BACKGROUND_COLOR = Ref("transparent")

const CONFIGHELP = """
### Extra Keyword Arguments
* `background_color = DEFAULT_BACKGROUND_COLOR[]`
* `scale::Int = 60.0` is the number of pixel per unit distance
* `xpad::Float64 = 2.5` is the padding space in x axis
* `ypad::Float64 = 1.5` is the padding space in y axis

##### axes
* `axes_text_color = DEFAULT_TEXT_COLOR[]`
* `axes_text_fontsize::Float64 = 16.0`
* `axes_num_of_xticks = 5`
* `axes_num_of_yticks = 5`
* `axes_x_offset::Float64 = 0.1`
* `axes_y_offset::Float64 = 0.06`
* `axes_unit::String = "μm"`

##### node
* `node_text_fontsize::Float64 = 16.0`
* `node_text_color = DEFAULT_TEXT_COLOR[]`
* `node_stroke_color = DEFAULT_LINE_COLOR[]`
* `node_stroke_linewidth = 1`
* `node_fill_color = DEFAULT_NODE_COLOR[]`

##### bond
* `bond_color = DEFAULT_LINE_COLOR[]`
* `bond_linewidth::Float64 = 1.0`

##### `blockade`
* `blockade_radius::Float64=0`, atoms within `blockade_radius` will be connected by edges.
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
* `grid_stroke_color="#AAAAAA"`
* `grid_stroke_width::Float64=1`
* `grid_stroke_style::String="dashed"`
"""

Base.@kwdef mutable struct AxisDisplayConfig
    axes_text_color::String = DEFAULT_TEXT_COLOR[]
    axes_text_fontsize::Float64 = 16.0
    axes_num_of_xticks::Int = 5
    axes_num_of_yticks::Int = 5
    axes_x_offset::Float64 = 0.1
    axes_y_offset::Float64 = 0.06
    axes_unit::String = "μm"
end

"""
    img_atoms(atoms::AtomList;
        colors = [DEFAULT_LINE_COLOR[], ...],
        texts = ["1", "2", ...],
        vectors = [],
        format = :svg,
        filename = nothing,
        kwargs...
        )
    img_atoms(bounded_lattice::BoundedLattice;
        colors = [DEFAULT_LINE_COLOR[], ...],
        texts = ["1", "2", ...],
        vectors = [],
        format = :svg,
        filename = nothing,
        kwargs...
        )
Plots `atoms` or `bounded_lattice.site_positions` with colors specified by `colors` and texts specified by `texts`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.
If you want to write this image to the disk without displaying it in a frontend, please try

```julia
julia> img_atoms(generate_sites(SquareLattice(), 5, 5); filename="test.png")
```

### Keyword Arguments
* `colors` is a vector of colors for nodes
* `texts` is a vector of string displayed on nodes
* `vectors` is a vector of (start_loc, end_loc) pair to specify a list of arrows.
* `format` can be `:svg`, `:pdf` or `:png`
* `filename` can be a filename string with suffix `.svg`, `.png` or `.pdf`

$CONFIGHELP
"""
function img_atoms(
        atoms::AtomList{2};
        # vertex features
        colors = nothing,
        texts = nothing,
        vectors = [],
        # IO
        format = :svg,
        filename = nothing,
        # config
        kwargs...
    )
    # config
    config = LatticeDisplayConfig(; kwargs...)

    #length(atoms) == 0 && return LuxorGraphPlot._draw(()->nothing, 100, 100; format, filename)
    edges = _edges(atoms, config.blockade_radius)
    atoms_config = config.blockade_radius > config.node_size ? extract_blockade_graph_display_config(config) : extract_atom_graph_display_config(config)
    diag = LuxorGraphPlot.diagram(map(a->a .* config.scale, atoms), edges; config=atoms_config)
    xpad, ypad = config.xpad * config.scale, config.ypad * config.scale
    with_nodes(diag; format, filename,
                    padding_bottom = ypad,
                    padding_left = xpad,
                    padding_right = xpad,
                    padding_top = ypad,
                    background = config.background_color
                ) do
        _viz_atoms(diag, colors, vectors .* config.scale, texts, config)
        xmin, xmax, ymin, ymax = minimum(a->a[1], atoms), maximum(a->a[1], atoms), minimum(a->a[2], atoms), maximum(a->a[2], atoms)
        _viz_axes(xmin-config.xpad/2, ymin-config.ypad/2, xmax+config.xpad/2, ymax+config.ypad/2, config.scale, axes_config)
    end
end
img_atoms(atoms::AtomList{1}; kwargs...) = img_atoms(padydim(atoms); kwargs...)
img_atoms(bounded_lattice::BoundedLattice{<:AbstractLattice{2},<:AbstractRegion{2}}; kwargs...) = img_atoms(bounded_lattice.site_positions; kwargs...)
img_atoms(bounded_lattice::BoundedLattice{<:AbstractLattice{1},<:AbstractRegion{1}}; kwargs...) = img_atoms(bounded_lattice.site_positions; kwargs...)


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
function _viz_axes(xmin, xmax, ymin, ymax, scale, config::AxisDisplayConfig)
    # the true coordinates
    xs = _LinRange(xmin, xmax, config.axes_num_of_xticks)
    ys = _LinRange(ymin, ymax, config.axes_num_of_yticks)
    coo(x, y) = Point((x, y)) * scale
    xys = [xs..., ys...]
    # the display coordinates of axes labels
    locs = [[coo(x, ymax+config.axes_y_offset) for x in xs]...,
        [coo(xmin-config.axes_x_offset, y) for y in ys]...]
    Luxor.@layer begin
        Luxor.sethue(config.axes_text_color)
        Luxor.fontsize(config.axes_text_fontsize)
        for (x, loc) in zip(xys, locs)
            Luxor.text("$(round(x; digits=2))$(config.axes_unit)", loc)
        end
    end
end

Base.@kwdef struct LatticeDisplayConfig
    # line, node and text
    background_color = DEFAULT_BACKGROUND_COLOR[]
    scale::Float64 = 60.0
    xpad::Float64 = 2.5
    ypad::Float64 = 1.5

    # axes
    axes_text_color = DEFAULT_LINE_COLOR[]  # NOTE: follow the line color!
    axes_text_fontsize::Float64 = 16.0
    axes_num_of_xticks::Int = 5
    axes_num_of_yticks::Int = 5
    axes_x_offset::Float64 = 0.1
    axes_y_offset::Float64 = 0.06
    axes_unit::String = "μm"

    # node
    node_size::Float64 = 0.25
    node_text_fontsize::Float64 = 16.0
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
function extract_atom_graph_display_config(config::LatticeDisplayConfig)
    GraphDisplayConfig(; background=config.background_color, fontsize=config.node_text_fontsize,
        vertex_shape=:circle, vertex_line_width=config.node_stroke_linewidth, vertex_line_style="solid",
        vertex_text_color=config.node_text_color, vertex_stroke_color=config.node_stroke_color, vertex_color=config.node_fill_color,
        vertex_size=config.node_size * config.scale, edge_color=config.bond_color, edge_line_width=config.bond_linewidth, edge_line_style="solid")
end
function extract_blockade_graph_display_config(config::LatticeDisplayConfig)
    GraphDisplayConfig(; background=config.background_color, fontsize=config.node_text_fontsize,
        vertex_shape=:circle, vertex_line_width=config.blockade_stroke_linewidth, vertex_line_style="dashed",
        vertex_text_color=config.node_text_color, vertex_stroke_color=config.blockade_stroke_color, vertex_color=config.blockade_fill_color,
        vertex_size=config.blockade_radius * config.scale, edge_color=config.bond_color, edge_line_width=config.bond_linewidth, edge_line_style="solid")
end

function _viz_atoms(diag, colors, vectors, texts, config::LatticeDisplayConfig)
    # show the blockade
    if config.blockade_style != "none"
        show_diagram(diag; config = extract_blockade_graph_display_config(config))
    end

    # show the arrows
    Luxor.@layer begin
        Luxor.sethue(config.arrow_color)
        Luxor.setline(config.arrow_linewidth)
        Luxor.setdash("solid")
        for v in vectors
            Luxor.arrow(Point(v[1])*config.unit, Point(v[2])*config.unit; arrowheadlength=config.arrow_head_length)
        end
    end

    # show the graph
    LuxorGraphPlot.show_diagram(diag; config=extract_atom_graph_display_config(config), texts,
        vertex_colors=colors,
        vertex_stroke_colors=nothing,
        vertex_text_colors=nothing,
        edge_colors=nothing
        )
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
        colors=[DEFAULT_LINE_COLOR[], ...],
        texts=["1", "2", ...],
        vectors=[],
        format=:svg,
        filename=nothing,
        kwargs...
        )

Draw a `maskedgrid` with colors specified by `colors` and texts specified by `texts`.
You will need a `VSCode`, `Pluto` notebook or `Jupyter` notebook to show the image.

### Keyword Arguments
* `colors` is a vector of colors for nodes
* `texts` is a vector of string displayed on nodes
* `vectors` is a vector of arrows
* `format` can be `:svg`, `:pdf` or `:png`
* `filename` can be a filename string with suffix `.svg`, `.png` or `.pdf`

$CONFIGHELP
"""
function img_maskedgrid(
        maskedgrid::MaskedGrid;
        format = :svg,
        filename = nothing,
        colors = nothing,
        texts = nothing,
        vectors = [],
        xpad = nothing,
        ypad = nothing,
        kwargs...,
    )
    atoms = padydim(collect_atoms(maskedgrid))
    xmin, ymin, xmax, ymax = LuxorGraphPlot.get_bounding_box(atoms)
    auto = config_plotting(atoms, xpad, ypad)

    config = LatticeDisplayConfig(; auto..., kwargs...)
    Dx, Dy = ((xmax-xmin)+2*auto.xpad)*config.scale*config.unit, ((ymax-ymin)+2*auto.ypad)*config.scale*config.unit
    transform(loc) = config.scale .* (loc[1]-xmin+auto.xpad, loc[2]-ymin+auto.ypad)
    LuxorGraphPlot._draw(Dx, Dy; format, filename) do
        LuxorGraphPlot.background(config.background_color)
        # show the grid
        _viz_grid(maskedgrid.xs, maskedgrid.ys; transform, unit=config.unit, auto.xpad, auto.ypad, xmin, ymin,
            xmax, ymax,
            color=config.grid_stroke_color, line_width=config.grid_stroke_width, line_style=config.grid_stroke_style)
        # show atoms
        length(atoms) > 0 && _viz_atoms(transform.(atoms), _edges(atoms, config.blockade_radius),
                colors, vectors, texts, config)
        # show the axes
        _viz_axes(; transform, xmin, ymin, xmax, ymax,
            axes_x_offset=config.axes_x_offset, axes_y_offset=config.axes_y_offset,
            axes_num_of_xticks=config.axes_num_of_xticks, axes_num_of_yticks=config.axes_num_of_yticks,
            axes_text_fontsize=config.axes_text_fontsize, axes_text_color=config.axes_text_color, axes_unit=config.axes_unit, unit=config.unit)
    end
end

function _viz_grid(xs, ys; transform, xpad, ypad, unit, xmin, ymin, xmax, ymax, color, line_width, line_style)
    _xmin, _xmax = (xmin - 0.3*xpad, xmax + 0.3*xpad)
    _ymin, _ymax = (ymin - 0.3*ypad, ymax + 0.3*ypad)
    coo(x, y) = Point(transform((x, y)))*unit
    for x in xs
        LuxorGraphPlot.draw_edge(coo(x, _ymin), coo(x, _ymax); color, line_width, line_style)
    end
    for y in ys
        LuxorGraphPlot.draw_edge(coo(_xmin, y), coo(_xmax, y); color, line_width, line_style)
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

        function Base.show(io::IO, m::$mime, bounded_lattice::BoundedLattice)
            Base.show(io, m, img_atoms(bounded_lattice.site_positions; format = $(QuoteNode(format))))
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

function Base.show(io::IO, mime::MIME"text/plain", atoms::AtomList)

    # convert list of tuples into two separate lists
    x = Float64[]
    y = Float64[]
    for coord in atoms
        push!(x, coord[1])
        
        # Could be the case that coordinate only has single x value
        # in which case we set the y value to 0
        if length(coord) == 1
            push!(y, 0)
        else
            push!(y, coord[2])
        end

    end

    plt = scatterplot(
        x, 
        y, 
        title = "Approximate Atom Positions",
        xlabel = "μm", 
        ylabel = "μm",
        marker = :circle)

    return show(io, mime, plt)
end