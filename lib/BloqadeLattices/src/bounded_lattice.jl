





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




function parallelepiped_region(lattice::AbstractLattice{D},M::Vararg{NTuple{D,Int},D};pbc::Bool=false) where D
    lat_vecs = lattice_vectors(lattice)
    T = eltype(lat_vecs[1])
    bounds =  zeros(T,D,D)
    
    for i in 1:D
        for j in 1:D
            bounds[:,i] .+= M[i][j] .* lat_vecs[j]
        end
    end
    region = Parallelepiped(bounds)
    
    return BoundedLattice(lattice,region,pbc)

end

dimension(lattice::BoundedLattice{L,C}) where {L,C} = dimension(lattice.lattice)
lattice_vectors(lattice::BoundedLattice{L,C}) where {L,C} = lattice_vectors(lattice.lattice)

# what behavior to produce when position is not found?
function get_position_index(pos,lattice::BoundedLattice{L,C}) where {L,C}
    j = searchsortedlast(Tuple(pos),lattice.site_positions)
    !(pos ≈ lattice.site_positions[j]) && error("$pos not contained in bounded lattice.")
    return j
end


