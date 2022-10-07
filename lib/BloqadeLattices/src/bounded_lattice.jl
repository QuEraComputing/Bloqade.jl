





# abstraction for a single tile of an infinite lattice

struct BoundedLattice{L<:AbstractLattice,R<:AbstractRegion} 
    lattice::L
    region::R
    site_positions::AtomList
    PBC::Bool
    function BoundedLattice(lattice::L,region::R,site_positions,PBC::Bool) where {L<:AbstractLattice,R<:AbstractRegion} 
        sort!(site_positions)

        return new{L,R}(lattice,region,AtomList(site_positions),PBC)
    end
end

function BoundedLattice(lattice::AbstractLattice{D},region::AbstractRegion{D},PBC::Bool=false) where D
    site_positions = generate_sites_in_region(lattice,region)
    return BoundedLattice(lattice,region,site_positions,PBC)
end




function parallelepiped_region(lattice::AbstractLattice{D},M::Vararg{NTuple{D,Int},D};PBC::Bool=false) where D
    lat_vecs = lattice_vectors(lattice)
    T = eltype(lat_vecs[1])
    bounds =  zeros(T,D,D)
    
    for i in 1:D
        for j in 1:D
            bounds[:,i] .+= M[i][j] .* lat_vecs[j]
        end
    end
    region = Parallelepiped(bounds)
    
    return BoundedLattice(lattice,region,PBC)

end

dimension(lattice::BoundedLattice{L,C}) where {L,C} = dimension(lattice.lattice)
lattice_vectors(lattice::BoundedLattice{L,C}) where {L,C} = lattice_vectors(lattice.lattice)

# what behavior to produce when position is not found?
function get_position_index(pos,lattice::BoundedLattice{L,C}) where {L,C}
    j = searchsortedlast(Tuple(pos),lattice.site_positions)
    !(pos ≈ lattice.site_positions[j]) && error("$pos not contained in bounded lattice.")
    return j
end


