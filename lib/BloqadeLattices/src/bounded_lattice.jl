





# abstraction for a single tile of an infinite lattice

struct BoundedLattice{L<:AbstractLattice,R<:AbstractRegion} 
    lattice::L
    region::R
    site_positions::AtomList
    pbc::Bool
    function BoundedLattice(lattice::L,region::R,site_positions,pbc::Bool) where {L<:AbstractLattice,R<:AbstractRegion} 
        sort!(site_positions)

        return new{L,R}(lattice,region,AtomList(site_positions),pbc)
    end
end

function BoundedLattice(lattice::AbstractLattice{D},region::AbstractRegion{D},pbc::Bool=false) where D
    site_positions = generate_sites_in_region(lattice,region)
    return BoundedLattice(lattice,region,site_positions,pbc)
end




function parallelepiped_region(lattice::AbstractLattice{D},M::Vararg{NTuple{D,Int},D};pbc::Bool=false,scale::Real=1) where D
    scale > 0 || error("scale must be a positive real value.")
    lat_vecs = lattice_vectors(lattice)
    lat_sites = lattice_sites(lattice)
    T = eltype(lat_vecs[1])
    scaled_bounds =  zeros(T,D,D)
    scaled_lat_vecs = Tuple(scale .* lat_vec for lat_vec in lat_vecs)
    scaled_lat_sites = Tuple(scale .* lat_site for lat_site in lat_sites)

    for i in 1:D
        for j in 1:D
            scaled_bounds[:,i] .+= M[i][j] .* scaled_lat_vecs[j]
        end
    end

    scaled_region = Parallelepiped(scaled_bounds)
    scaled_lattice = GeneralLattice(scaled_lat_vecs,scaled_lat_sites)
    
    return BoundedLattice(scaled_lattice,scaled_region,pbc)

end

dimension(lattice::BoundedLattice{L,C}) where {L,C} = dimension(lattice.lattice)
lattice_vectors(lattice::BoundedLattice{L,C}) where {L,C} = lattice_vectors(lattice.lattice)
Base.length(lattice::BoundedLattice{L,C}) where {L,C} = length(lattice.site_positions)

# what behavior to produce when position is not found?
function get_position_index(pos,lattice::BoundedLattice{L,C}) where {L,C}
    j = searchsortedlast(Tuple(pos),lattice.site_positions)
    !(pos ≈ lattice.site_positions[j]) && error("$pos not contained in bounded lattice.")
    return j
end

distance(lat::BoundedLattice,x,y) = lat.pbc ? distance(lat.region,x,y) : distance(x,y)
