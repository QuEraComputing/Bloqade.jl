struct CuSparseDeviceMatrixCSR{Tv} <: AbstractCuSparseMatrix{Tv}
    rowPtr::CuDeviceVector{Cint, AS.Global}
    colVal::CuDeviceVector{Cint, AS.Global}
    nzVal::CuDeviceVector{Tv, AS.Global}
    dims::NTuple{2, Int}
    nnz::Cint
end

Base.size(H::CuSparseDeviceMatrixCSR) = H.dims

function Adapt.adapt_structure(to::CUDA.Adaptor, x::CuSparseMatrixCSR{Tv}) where Tv
    CuSparseDeviceMatrixCSR(cudaconvert(x.rowPtr), cudaconvert(x.colVal), cudaconvert(x.nzVal), x.dims, x.nnz)
end

function Base.show(io::IO, ::MIME"text/plain", A::CuSparseDeviceMatrixCSR)
    println(io, "$(length(A))-element device sparse matrix CSR at:")
    println(io, "  rowPtr $(pointer(A.rowPtr))")
    println(io, "  colVal $(pointer(A.colVal))")
    print(io, "  nzVal $(pointer(A.nzVal))")
end

function Adapt.adapt_structure(to::CUDA.Adaptor, t::XTerm)
    XTerm(t.nsites, cudaconvert(t.Ωs), cudaconvert(t.ϕs))
end

function Adapt.adapt_structure(to::CUDA.Adaptor, t::ZTerm)
    ZTerm(t.nsites, cudaconvert(t.Δs))
end

function Adapt.adapt_structure(to::CUDA.Adaptor, t::RydInteract)
    RydInteract(cudaconvert(t.C), cudaconvert(t.atoms))
end

function Adapt.adapt_structure(to::CUDA.Adaptor, t::Hamiltonian)
    Hamiltonian(cudaconvert.(t.terms))
end

function Adapt.adapt_structure(to, r::RydbergReg{N}) where {N}
    return RydbergReg{N}(adapt(to, r.state), adapt(to, r.subspace))
end

Adapt.adapt_structure(to, s::Subspace) = Subspace(s.map, adapt(to, s.subspace_v))
