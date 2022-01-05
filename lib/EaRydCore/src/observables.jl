module Op

using Yao
const n = Yao.ConstGate.P1

end

function rydberg_density(reg::AbstractRegister)
    nsites = nqubits(reg)
    return [expect(put(nsites, i=>Op.n), reg) for i in 1:nsites]    
end

function correlations(reg::AbstractRegister, bonds::Vector{<:Tuple})
    return map(bonds) do bond
        cor = put(nsites, bond=>repeat(length(bond), Op.n))
        expect(cor, info.reg)
    end
end
