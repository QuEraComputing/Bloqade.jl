# bond_spins.jl
#
# Defines the spatial lattice; supports 1D and 2D square, open and periodic boundaries
# The main data structure is the bond-spin index array, bond_spin[nBond,2]

lattice_bond_spins((nX,)::Tuple{Int}, pbc::Bool=true) = lattice_bond_spins(nX, pbc)

function lattice_bond_spins(nX::Int, pbc::Bool=true)
    Nb = pbc ? nX : nX-1
    Ns = nX

    bond_spin = [(0, 0) for _ in 1:Nb]
    for i in 1:Nb
        bond_spin[i] = (i, i+1)
    end

    if pbc
        bond_spin[end] = (Nb, 1)
    end

    return bond_spin, Ns, Nb
end

# 2D square lattice
function lattice_bond_spins(nX::NTuple{2, Int}, pbc::Bool=true)
    (nX[1] == nX[2]) || throw(ArgumentError("Currently only supports square lattices!"))
    L = nX[1]

    Ns = L * L
    Nb = pbc ? 2 * Ns : 2 * (Ns - L)
    bond_spin = [(0, 0) for _ in 1:Nb]

    if pbc
        @inbounds for i in 1:Ns
            # horizontal
            bond1 = 2 * i - 1

            if mod(i, L) != 0
                bond_spin[bond1] = (i, i + 1)
            else
                bond_spin[bond1] = (i, i + 1 - L)
            end

            # vertical
            bond2 = 2 * i

            if i <= (Ns - L)
                bond_spin[bond2] = (i, i + L)
            else
                bond_spin[bond2] = (i, i + L - Ns)
            end
        end
    else
        # horizontal
        cnt = 1
        @inbounds for i in 1:div(Nb, 2)
            if mod(i, L - 1) != 0
                bond_spin[i] = (cnt, cnt + 1)
                cnt += 1
            else
                bond_spin[i] = (cnt, cnt + 1)
                cnt += 2
            end
        end

        # vertical
        cnt = 1
        @inbounds for i in (div(Nb, 2)+1):Nb
            bond_spin[i] = (cnt, cnt + L)
            cnt += 1
        end
    end

    return bond_spin, Ns, Nb
end
