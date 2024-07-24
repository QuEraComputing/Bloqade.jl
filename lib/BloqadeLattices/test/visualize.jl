using Test
using BloqadeLattices
using LuxorGraphPlot
using Documenter

@testset "visualize (image)" begin
    BloqadeLattices.darktheme!()
    lt = generate_sites(KagomeLattice(), 5, 5, scale = 1.5)
    blt = parallelepiped_region(KagomeLattice(), (5,0),(0,5);scale=1.5)
    grd = make_grid(lt[2:end-1])
    unitvectors(lattice::AbstractLattice, scale::Real) = [((0.0, 0.0), v .* scale) for v in lattice_vectors(lattice)]
    @test img_atoms(lt; vectors = unitvectors(KagomeLattice(), 1.5)) isa LuxorGraphPlot.Drawing
    # different colors
    @test img_atoms(lt; colors = nothing) isa LuxorGraphPlot.Drawing
    @test img_atoms(lt; node_fill_color = "red") isa LuxorGraphPlot.Drawing
    @test img_atoms(lt; colors = fill("blue", length(lt))) isa LuxorGraphPlot.Drawing
    @test img_atoms(lt; colors = ByDensity(rand(0:10, length(lt)); vmax = 10)) isa LuxorGraphPlot.Drawing
    @test img_atoms(blt; colors = nothing) isa LuxorGraphPlot.Drawing
    @test img_atoms(blt; node_fill_color = "red") isa LuxorGraphPlot.Drawing
    @test img_atoms(blt; colors = fill("blue", length(lt))) isa LuxorGraphPlot.Drawing
    @test img_atoms(blt; colors = ByDensity(rand(0:10, length(lt)); vmax = 10)) isa LuxorGraphPlot.Drawing
    @test img_maskedgrid(grd) isa LuxorGraphPlot.Drawing
    @test show(IOBuffer(), MIME"image/svg+xml"(), grd) === nothing
    @test show(IOBuffer(), MIME"image/svg+xml"(), lt) === nothing
    @test show(IOBuffer(), MIME"image/png"(), grd) === nothing
    @test show(IOBuffer(), MIME"image/png"(), lt) === nothing

    BloqadeLattices.lighttheme!()
    @test img_atoms(lt; colors = nothing) isa LuxorGraphPlot.Drawing
    @test show(IOBuffer(), MIME"image/svg+xml"(), lt) === nothing
end


@testset "visualize (text)" begin

    @testset "single x value" begin
        lt = generate_sites(ChainLattice(), 10, scale=5.74)
        show(stdout, MIME"text/plain"(), lt)
    end 

    @testset "x and y values" begin
        lt = generate_sites(SquareLattice(), 5, 5, scale=5.74)
        show(stdout, MIME"text/plain"(), lt)
    end

end