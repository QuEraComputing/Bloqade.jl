# D is the dimensionality
abstract type AbstractLattice{D} end

# K is the sites in a unit cell,
# T is normally floating point numbers for locations,
struct BravaisLattice{D, K, T<:Real}
    bases::NTuple{D,NTuple{D,T}}
    sites::NTuple{K, NTuple{D,T}}
end

function generate_sites(bl::BravaisLattice{D,K,T}, repeats...) where {D,K,T}
    @assert D > 0
    @assert all(>=(0), repeats)
    @assert length(repeats) == D
    locations = NTuple{D,T}[]  # we might want to avoid using `push!` later.
    for ci in CartesianIndices(repeats)
        baseloc = mapreduce(i->(ci.I[i]-1) .* bl.bases[i], (x, y) -> x .+ y, 1:D)
        for siteloc in bl.sites
            push!(locations, baseloc .+ siteloc)
        end
    end
    return locations
end

struct GeneralLattice{D,K,T} <: AbstractLattice{D}
    vectors::NTuple{D,NTuple{D,T}}
    sites::NTuple{K,NTuple{D,T}}
end
GeneralLattice(vectors, sites) = GeneralLattice((vectors...,), (sites...,))
lattice_vectors(gl::GeneralLattice) = gl.vectors
lattice_sites(gl::GeneralLattice) = gl.sites

struct HoneycombLattice <: AbstractLattice{2} end
lattice_vectors(::HoneycombLattice) = ((1.0, 0.0), (0.5, 0.5*sqrt(3)))
lattice_sites(::HoneycombLattice) = ((0.0, 0.0), (0.5, 0.5/sqrt(3)))

struct SquareLattice <: AbstractLattice{2} end
lattice_vectors(::SquareLattice) = ((1.0, 0.0), (0.0, 1.0))
lattice_sites(::SquareLattice) = ((0.0, 0.0),)

struct TriangularLattice <: AbstractLattice{2} end
lattice_vectors(::TriangularLattice) = ((1.0, 0.0), (0.5, 0.5*sqrt(3)))
lattice_sites(::TriangularLattice) = ((0.0, 0.0),)

struct ChainLattice <: AbstractLattice{1} end
lattice_vectors(::ChainLattice) = ((1.0,),)
lattice_sites(::ChainLattice) = ((0.0,),)

struct LiebLattice <: AbstractLattice{2} end
lattice_vectors(::LiebLattice) = ((1.0, 0.0), (0.0, 1.0))
lattice_sites(::LiebLattice) = ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5))

struct KagomeLattice <: AbstractLattice{2} end
lattice_vectors(::KagomeLattice) = ((1.0, 0.0), (0.5, 0.5*sqrt(3)))
lattice_sites(::KagomeLattice) = ((0.0, 0.0), (0.25, 0.25*sqrt(3)), (0.75, 0.25*sqrt(3)))

# lattice -> bravais lattice
bravais(lt::AbstractLattice) = BravaisLattice((lattice_vectors(lt)...,), (lattice_sites(lt)...,))
generate_sites(lt::AbstractLattice, nrepeats...) = generate_sites(bravais(lt), nrepeats...)

############ manipulate sites ###############
# offset sites
function offset_axes(sites::AbstractVector{NTuple{D, T}}, offsets...) where {D, T}
    @assert length(offsets) == D
    return map(x->ntuple(i->x[i]+offsets[i], D), sites)
end

# dropout sites
function random_dropout(sites::AbstractVector{NTuple{D, T}}, prob::Real) where {D, T}
    return sites[rand(length(sites)) .> prob]
end

# filter out sites out of bounds
function clip_axes(sites::AbstractVector{NTuple{D, T}}, bounds...) where {D, T}
    @assert length(bounds) == D
    @assert all(x->length(x) == 2, bounds)
    return filter(x->all(i->bounds[i][1] <= x[i] <= bounds[i][2], 1:D), sites)
end
clip_axes(args...) = ls -> clip_axes(ls, args...)
offset_axes(args...) = ls -> offset_axes(ls, args...)
random_dropout(prob::Real) = ls -> random_dropout(ls, prob)

############ manipulate grid ###############
struct MaskedGrid{T}
    xs::Vector{T}
    ys::Vector{T}
    mask::Matrix{Bool}
end

# create `MaskedGrid` from the locations.
function make_grid(sites::AbstractVector{NTuple{1, T}}; atol=10*eps(T)) where {T}
    make_grid(padydim.(sites); atol=atol)
end
padydim(x::Tuple{T}) where T = (x[1], zero(T))
function make_grid(sites::AbstractVector{NTuple{2, T}}; atol=10*eps(T)) where {T}
    xs = sort!(approximate_unique(getindex.(sites, 1), atol))
    ys = sort!(approximate_unique(getindex.(sites, 2), atol))
    ixs = map(s->findfirst(==(s[1]), xs), sites)
    iys = map(s->findfirst(==(s[2]), ys), sites)
    m, n = length(xs), length(ys)
    mask = zeros(Bool, m, n)
    for (ix, iy) in zip(ixs, iys)
        mask[ix, iy] = true
    end
    return MaskedGrid(xs, ys, mask)
end

# return `(uxs, ixs)``, where `uxs` is the unique x-coordinates, `ixs` the mapping from the index in `xs` to the index in `uxs`.
function approximate_unique(xs::AbstractVector{T}, atol) where T
    uxs = T[]
    for x in xs
        found = false
        for ux in uxs
            if isapprox(x, ux; atol=atol)
                found = true
                break
            end
        end
        if !found
            push!(uxs, x)
        end
    end
    return uxs
end

# get locations in order
function locations(mg::MaskedGrid)
    map(ci->(mg.xs[ci.I[1]], mg.ys[ci.I[2]]), findall(mg.mask))
end

# TODO
# pseudo-lattices,
# image/svg output (maybe),