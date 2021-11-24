export AbstractLattice, BravaisLattice, HoneycombLattice, SquareLattice, TriangularLattice, ChainLattice, LiebLattice, KagomeLattice
export bravais, generate_sites, offsetaxes, clipaxes, latticesites, latticevectors
export MaskedGrid, makegrid, locations

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

struct HoneycombLattice <: AbstractLattice{2} end
latticevectors(::HoneycombLattice) = ((1.0, 0.0), (0.5, 0.5*sqrt(3)))
latticesites(::HoneycombLattice) = ((0.0, 0.0), (0.5, 0.5/sqrt(3)))

struct SquareLattice <: AbstractLattice{2} end
latticevectors(::SquareLattice) = ((1.0, 0.0), (0.0, 1.0))
latticesites(::SquareLattice) = ((0.0, 0.0),)

struct TriangularLattice <: AbstractLattice{2} end
latticevectors(::TriangularLattice) = ((1.0, 0.0), (0.5, 0.5*sqrt(3)))
latticesites(::TriangularLattice) = ((0.0, 0.0),)

struct ChainLattice <: AbstractLattice{1} end
latticevectors(::ChainLattice) = ((1.0,),)
latticesites(::ChainLattice) = ((0.0,),)

struct LiebLattice <: AbstractLattice{2} end
latticevectors(::LiebLattice) = ((1.0, 0.0), (0.0, 1.0))
latticesites(::LiebLattice) = ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5))

struct KagomeLattice <: AbstractLattice{2} end
latticevectors(::KagomeLattice) = ((1.0, 0.0), (0.5, 0.5*sqrt(3)))
latticesites(::KagomeLattice) = ((0.0, 0.0), (0.25, 0.25*sqrt(3)), (0.75, 0.25*sqrt(3)))

# lattice -> bravais lattice
bravais(lt::AbstractLattice) = BravaisLattice((latticevectors(lt)...,), (latticesites(lt)...,))
generate_sites(lt::AbstractLattice, nrepeats...) = generate_sites(bravais(lt), nrepeats...)

############ manipulate sites ###############
# offset sites
function offsetaxes(sites::AbstractVector{NTuple{D, T}}, offsets...) where {D, T}
    @assert length(offsets) == D
    return map(x->ntuple(i->x[i]+offsets[i], D), sites)
end

# filter out sites out of bounds
function clipaxes(sites::AbstractVector{NTuple{D, T}}, bounds...) where {D, T}
    @assert length(bounds) == D
    @assert all(x->length(x) == 2, bounds)
    return filter(x->all(i->bounds[i][1] <= x[i] <= bounds[i][2], 1:D), sites)
end

############ manipulate grid ###############
struct MaskedGrid{T}
    xs::Vector{T}
    ys::Vector{T}
    mask::Matrix{Bool}
end

# create `MaskedGrid` from the locations.
function makegrid(sites::AbstractVector{NTuple{2, T}}; atol=10*eps(T)) where {T}
    xs, ixs = approximate_unique(getindex.(sites, 1), atol)
    ys, iys = approximate_unique(getindex.(sites, 2), atol)
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
    ixs = Vector{Int}(undef, length(xs))
    for (i, x) in enumerate(xs)
        found = false
        for (k, ux) in enumerate(uxs)
            if isapprox(x, ux; atol=atol)
                ixs[i] = k
                found = true
                break
            end
        end
        if !found
            push!(uxs, x)
            ixs[i] = length(uxs)
        end
    end
    return uxs, ixs
end

# get locations in order
function locations(mg::MaskedGrid)
    map(ci->(mg.xs[ci.I[1]], mg.ys[ci.I[2]]), findall(mg.mask))
end

# TODO
# pseudo-lattices,
# image/svg output,
# use KDTree.