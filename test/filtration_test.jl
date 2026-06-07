@testset "Filtration" begin
    # Basic filtration
    s1 = Simplex(1)
    s2 = Simplex(2)
    s3 = Simplex(1, 2)
    f = Filtration([s1, s2, s3], [0.0, 0.0, 1.5])

    @test length(f) == 3
    @test filtration_value(f, 1) == 0.0
    @test filtration_value(f, 3) == 1.5

    # simplices_below
    below = simplices_below(f, 0.0)
    @test length(below) == 2  # two vertices

    below2 = simplices_below(f, 2.0)
    @test length(below2) == 3

    # complex_at
    sc = complex_at(f, 0.5)
    @test has_simplex(sc, Simplex(1))
    @test has_simplex(sc, Simplex(2))
    @test !has_simplex(sc, Simplex(1, 2))

    # Vietoris-Rips
    points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
    vrf = VietorisRipsFiltration(points, 2, 2.0)
    @test length(vrf) >= 6  # 3 vertices + 3 edges + possibly 1 triangle

    # VR with matrix input
    pts_matrix = Float64[0.0 0.0; 1.0 0.0; 0.0 1.0]
    vrf2 = VietorisRipsFiltration(pts_matrix, 1, 2.0)
    @test length(vrf2) >= 5  # 3 vertices + 3 edges

    # Empty input
    f_empty = Filtration(Simplex[], Float64[])
    @test length(f_empty) == 0

    # euclidean_dist
    @test euclidean_dist([0.0, 0.0], [3.0, 4.0]) ≈ 5.0

    # Collinear points (degenerate)
    col_points = [[0.0], [1.0], [2.0], [3.0]]
    vrf_col = VietorisRipsFiltration(col_points, 2, 2.0)
    @test length(vrf_col) >= 4  # at least 4 vertices
end
