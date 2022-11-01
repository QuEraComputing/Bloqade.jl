
function two_body_interaction_matrix(f, atoms)

    natoms = length(atoms)
    mat = zeros(natoms,natoms)

    for i in 1:natoms
        for j in i+1:natoms
            @inbounds mat[i,j] = f(atoms[i],atoms[j])
        end
    end

    return mat
end 

function rydberg_interaction_matrix(atoms, C::Real)
    return two_body_interaction_matrix(atoms) do x,y
        return C/distance(x,y)^6
    end
end


### implementation for bounded lattices


function rydberg_interaction_matrix(lat::BoundedLattice{L,R},C::Real) where {L,R}
    if lat.PBC
        return two_body_interaction_matrix(lat.site_positions) do x,y
            return C/distance(lat.region,x,y)^6
        end 
    else
        return two_body_interaction_matrix(lat.site_positions) do x,y
            return C/distance(x,y)^6
        end
    end
end

#### TODO: Add generator versions of this code. return (i,j,Vij)