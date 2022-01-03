using CSV
using DataFrames
using DelimitedFiles
using EaRydCore
using EaRydODE

function save_state(r::RydbergReg, L, graph_index, T_index, R)
    df = DataFrame(:state=>vec(r.state), :subspace=>r.subspace.subspace_v)
    CSV.write("state-L-$L-graph-$graph_index-T-$T_index-R-$R.csv", df)
end

function save_ratio(ratio::Real, added_ratio::Real, L, graph_index, T_index, R)
    writedlm("ratio-L-$L-graph-$graph_index-T-$T_index-R-$R.txt", [ratio, added_ratio])
end

function read_lattice(L::Int, graph_index::Int)
    graph_data = readdlm(joinpath(@__DIR__, "graph/L$L.dat"))
    lattice_cfg = graph_data[graph_index, :][4:end]
    len = Int(sqrt(length(lattice_cfg)))
    return reshape(lattice_cfg, len, len)
end

function lattice_to_atoms(lattice::Matrix)
    return map(findall(isequal(1), lattice)) do x
        RydAtom(Tuple(x))
    end
end