struct PrecisionAdaptor{P} end

PrecisionAdaptor(P::Type) = PrecisionAdaptor{P}()
Adapt.adapt_storage(::PrecisionAdaptor{P}, x::Real) where P = P(x)
Adapt.adapt_storage(::PrecisionAdaptor{P}, x::Complex) where P = Complex{P}(x)
Adapt.adapt_storage(::PrecisionAdaptor{P}, x::Array) where P = convert(Array{P}, x)
Adapt.adapt_storage(::PrecisionAdaptor{P}, x::Array{<:Complex}) where P = convert(Array{Complex{P}}, x)

function storage_size(cache::DiscreteEmulationCache)
    return storage_size(cache.H)
end

function storage_size(H::SparseMatrixCSC)
    return sizeof(H.colptr) + sizeof(H.rowval) + sizeof(H.nzval)
end

storage_size(r::RydbergReg) = sizeof(r.state) + sizeof(r.subspace.subspace_v) + sizeof(r.subspace.map)
storage_size(r::Yao.ArrayReg) = sizeof(r.state)
storage_size(H::Array) = sizeof(H)
storage_size(x) = sizeof(x) # fallback to sizeof

get_space(r::Yao.ArrayReg) = fullspace
get_space(r::RydbergReg) = r.subspace

num_zero_term(t::Hamiltonian) = count(iszero, t.terms)
num_zero_term(t::AbstractTerm) = iszero(t) ? 1 : 0
