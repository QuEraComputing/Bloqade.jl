using Test
using BloqadeExpr, YaoBlocks, BitBasis
using SparseArrays

@testset "getindex" begin
    pb = rydberg_h([(5*randn(), 5*randn()) for i=1:5]; Ω=0.3, Δ=0.5)
    mpb = mat(pb)
    allpass = true
    for i=basis(pb), j=basis(pb)
        allpass &= pb[i, j] == mpb[Int(i)+1, Int(j)+1]
    end
    @test allpass

    allpass = true
    for j=basis(pb)
        allpass &= vec(pb[:, j]) == mpb[:, Int(j)+1]
        allpass &= vec(pb[j,:]) == mpb[Int(j)+1,:]
        allpass &= vec(pb[:, EntryTable([j], [1.0+0im])]) == mpb[:, Int(j)+1]
        allpass &= vec(pb[EntryTable([j], [1.0+0im]),:]) == mpb[Int(j)+1,:]
    end
    @test allpass
end