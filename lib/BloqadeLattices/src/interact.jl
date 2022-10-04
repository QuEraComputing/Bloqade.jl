



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


function rydberg_interaction_matrix(lat::BoundedLattice,C::Real)
    if lat.PBC
        return two_body_interaction_matrix(lat.lattice_sites,(x,y)->C/distance(lat.region,x,y)^6)
    else
        return two_body_interaction_matrix(lat.lattice_sites,(x,y)->C/distance(x,y)^6)
    end
end

#### TODO: Add generator versions of this code. return (i,j,Vij)