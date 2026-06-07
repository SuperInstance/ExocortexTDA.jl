@testset "Betti" begin
    # Betti numbers from barcode at threshold
    bc = PersistenceBarcode([(0.0, 1.0), (0.0, 2.0), (0.5, Inf)], [0, 0, 1])
    bn = betti_numbers(bc, 0.5)
    @test get(bn, 0, 0) == 2  # both H_0 bars alive at 0.5
    @test get(bn, 1, 0) == 1  # H_1 bar alive (infinite)

    # At threshold 1.5
    bn2 = betti_numbers(bc, 1.5)
    @test get(bn2, 0, 0) == 1  # only second bar alive

    # Betti numbers from filtration
    points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
    vrf = VietorisRipsFiltration(points, 2, 5.0)
    bn_f = betti_numbers(vrf, 0.01)
    @test get(bn_f, 0, 0) == 3  # three connected components at near-zero

    # Betti curve
    thresholds = [0.0, 0.5, 1.0, 1.5, 2.0]
    curve = betti_curve(bc, thresholds)
    @test length(curve[0]) == 5
    @test curve[0][1] == 2  # at 0.0, two H_0 bars alive

    # Betti curve from filtration
    curve_f = betti_curve(vrf, thresholds)
    @test haskey(curve_f, 0)

    # Simple complex: one triangle → β₀=1, β₁=0
    sc = SimplicialComplex()
    add_simplex!(sc, Simplex([1, 2, 3]))
    bn_sc = betti_numbers(sc)
    @test get(bn_sc, 0, 0) == 1

    # Hollow triangle (no 2-simplex): β₀=1, β₁=1
    sc2 = SimplicialComplex()
    add_simplex!(sc2, Simplex(1, 2))
    add_simplex!(sc2, Simplex(2, 3))
    add_simplex!(sc2, Simplex(1, 3))
    # Make sure we don't add the triangle
    bn_sc2 = betti_numbers(sc2)
    @test get(bn_sc2, 0, 0) == 1
end
