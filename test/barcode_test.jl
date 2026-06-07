using ExocortexTDA
using Test

@testset "Barcode" begin
    # Basic barcode
    bc = PersistenceBarcode([(0.0, 1.0), (0.0, 2.0), (0.5, Inf)], [0, 0, 1])
    @test num_bars(bc) == 3
    @test num_bars(bc, 0) == 2
    @test num_bars(bc, 1) == 1

    # Bars in dimension
    @test length(bars_in_dim(bc, 0)) == 2
    @test length(bars_in_dim(bc, 1)) == 1

    # Persistence diagram from barcode
    pd = PersistenceDiagram(bc)
    @test num_points(pd) == 3

    # Compute persistence from filtration
    points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
    vrf = VietorisRipsFiltration(points, 2, 5.0)
    bc2 = compute_persistence(vrf)
    @test num_bars(bc2) > 0
    # Should have H_0 bars
    @test num_bars(bc2, 0) >= 1

    # ASCII plot doesn't error
    ascii_str = plot_ascii(bc)
    @test occursin("Dimension", ascii_str)

    # ASCII diagram plot
    diag_str = plot_ascii(pd)
    @test typeof(diag_str) == String
    @test length(diag_str) > 0

    # Empty barcode
    bc_empty = PersistenceBarcode(Tuple{Float64,Float64}[], Int[])
    @test num_bars(bc_empty) == 0
    @test occursin("empty", plot_ascii(bc_empty))

    # Compute persistence from SimplicialComplex
    sc = SimplicialComplex()
    add_simplex!(sc, Simplex([1, 2, 3]))
    bc3 = compute_persistence(sc)
    @test num_bars(bc3) >= 1
end
