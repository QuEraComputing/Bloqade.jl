



function two_body_interaction_matrix(atoms,f)
    natoms = length(atoms)
    mat = zeros(natoms,natoms)

    for i in 1:natoms
        for j in i+1:natoms
            @inbounds mat[i,j] = f(atoms[i],atoms[j])
        end
    end

    return mat
end 

function rydberg_interaction_matrix(atoms,C::Real)
    return two_body_interaction_matrix(atoms,(x,y)->C/distance(x,y)^6)
end


### implementation for bounded lattices


function rydberg_interaction_matrix(lat::BoundedLattice{L,R},C::Real) where {L,R}
    return two_body_interaction_matrix(lat.site_positions,(x,y)->C/distance(lat,x,y)^6)
end

#### TODO: Add generator versions of this code. return (i,j,Vij)