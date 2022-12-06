"""
    two_body_interaction_matrix(f, atoms)

Generates an interaction matrix given a function `f` that accepts two atom positions at a time 
and an indexable iterable `atoms` containing atom positions.

See also [`rydberg_interaction_matrix`](@ref)

```jldoctest; setup=:(using BloqadeLattices)
julia> atoms = [(0.0,), (1.0,), (2.0,), (3.0,)]; # 1D chain, can also be AtomList

julia> two_body_interaction_matrix(atoms) do x,y return 1/distance(x,y) end
4×4 UpperTriangular{Float64, Matrix{Float64}}:
 0.0  1.0  0.5  0.333333
  ⋅   0.0  1.0  0.5
  ⋅    ⋅   0.0  1.0
  ⋅    ⋅    ⋅   0.0
```
"""
function two_body_interaction_matrix(f, atoms)

    natoms = length(atoms)
    mat = zeros(natoms,natoms)

    for i in 1:natoms
        for j in i+1:natoms
            @inbounds mat[i,j] = f(atoms[i],atoms[j])
        end
    end

    return UpperTriangular(mat)
end 

"""
    rydberg_interaction_matrix(atoms, C::Real)
    rydberg_interaction_matrix(lattice::BoundedLattice{L,R},C::Real)

Generate the interaction matrix given an indexable iterable `atoms` containg atom positions and the
Rydberg interaction constant `C`.

A `BoundedLattice` can be used in place of `atoms` which generates the Rydberg interaction matrix
for the lattice, factoring in Periodic Boundary Conditions.

See also [`two_body_interaction_matrix`](@ref)

```jldoctest; setup=:(using BloqadeLattices)
julia> atoms = [(0.0,), (1.0,), (2.0,), (3.0,)]; # 1D chain of atoms

julia> rydberg_interaction_matrix(atoms, 2π * 862690) # provide Rydberg constant
4×4 UpperTriangular{Float64, Matrix{Float64}}:
 0.0  5.42044e6  84694.4         7435.45
  ⋅   0.0            5.42044e6  84694.4
  ⋅    ⋅             0.0            5.42044e6
  ⋅    ⋅              ⋅             0.0

julia> bl = parallelepiped_region(SquareLattice(),(2,0),(0,2);pbc=true); 

julia> rydberg_interaction_matrix(bl, 2π * 862690)
4×4 UpperTriangular{Float64, Matrix{Float64}}:
 0.0  5.42044e6  5.42044e6  6.77555e5
  ⋅   0.0        6.77555e5  5.42044e6
  ⋅    ⋅         0.0        5.42044e6
  ⋅    ⋅          ⋅         0.0
```
"""
function rydberg_interaction_matrix(atoms, C::Real)
    return two_body_interaction_matrix(atoms) do x,y
        return C/distance(x,y)^6
    end
end


### implementation for bounded lattices


function rydberg_interaction_matrix(lattice::BoundedLattice{L,R},C::Real) where {L,R}
    if lattice.pbc
        return two_body_interaction_matrix(lattice.site_positions) do x,y
            return C/distance(lattice.region,x,y)^6
        end 
    else
        return two_body_interaction_matrix(lattice.site_positions) do x,y
            return C/distance(x,y)^6
        end
    end
end

#### TODO: Add generator versions of this code. return (i,j,Vij)