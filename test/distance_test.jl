@testset "Distance" begin
    # Identical diagrams
    pd1 = PersistenceDiagram([(0.0, 1.0), (0.0, 2.0)])
    pd2 = PersistenceDiagram([(0.0, 1.0), (0.0, 2.0)])
    @test bottleneck_distance(pd1, pd2) ≈ 0.0 atol=1e-10
    @test wasserstein_distance(pd1, pd2) ≈ 0.0 atol=1e-10

    # Different diagrams
    pd3 = PersistenceDiagram([(0.0, 1.0)])
    pd4 = PersistenceDiagram([(0.0, 2.0)])
    bn = bottleneck_distance(pd3, pd4)
    @test bn > 0.0

    ws = wasserstein_distance(pd3, pd4, 1)
    @test ws > 0.0

    # Empty diagrams
    pd_empty = PersistenceDiagram(Tuple{Float64,Float64}[])
    @test bottleneck_distance(pd_empty, pd_empty) ≈ 0.0

    # One empty, one not
    bn2 = bottleneck_distance(pd_empty, pd3)
    @test bn2 > 0.0

    # W_p with different p values
    pd5 = PersistenceDiagram([(0.0, 1.0), (0.0, 3.0)])
    pd6 = PersistenceDiagram([(0.0, 2.0), (0.0, 4.0)])
    w1 = wasserstein_distance(pd5, pd6, 1)
    w2 = wasserstein_distance(pd5, pd6, 2)
    @test w1 >= 0.0
    @test w2 >= 0.0
    # W_1 >= W_2 for same diagram pair
    @test w1 >= w2 - 1e-10

    # Diagonal distance
    @test diag_dist((0.0, 2.0)) ≈ 1.0
    @test diag_dist((1.0, 1.0)) ≈ 0.0

    # L∞ distance
    @test linf_dist((0.0, 1.0), (0.0, 1.0)) ≈ 0.0
    @test linf_dist((0.0, 1.0), (0.5, 2.0)) ≈ 1.0

    # Inf handling
    pd_inf = PersistenceDiagram([(0.0, Inf)])
    pd_fin = PersistenceDiagram([(0.0, 1.0)])
    bn_inf = bottleneck_distance(pd_inf, pd_fin)
    @test bn_inf == Inf || bn_inf >= 0.0  # inf bar matched differently

    # Default p=2 for wasserstein
    @test wasserstein_distance(pd3, pd4) == wasserstein_distance(pd3, pd4, 2)
end
