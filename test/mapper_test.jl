using ExocortexTDA
using Test

@testset "Mapper" begin
    # Simple mapper with custom filter
    points = [[0.0], [0.1], [0.9], [1.0], [2.0], [2.1]]

    # Identity filter
    mg = mapper_graph(points, p -> p[1], 3, 0.3, default_cluster)
    @test num_nodes(mg) >= 1
    @test num_edges(mg) >= 0
    @test typeof(mg) == MapperGraph

    # Symbol filter: :distance_to_mean
    mg2 = mapper_graph(points, :distance_to_mean, 2, 0.3, default_cluster)
    @test num_nodes(mg2) >= 1

    # Symbol filter: :eccentricity
    mg3 = mapper_graph(points, :eccentricity, 2, 0.3, default_cluster)
    @test num_nodes(mg3) >= 1

    # Matrix input (proper 2D matrix)
    pts_mat = Float64[0.0 0.0; 0.1 0.1; 0.9 0.9; 1.0 1.0; 2.0 2.0; 2.1 2.1]
    mg4 = mapper_graph(pts_mat, p -> p[1], 2, 0.3, default_cluster)
    @test num_nodes(mg4) >= 1

    # Default cluster function
    clusters = default_cluster(points, [1, 2, 3, 4])
    @test length(clusters) >= 1
    @test all(c -> !isempty(c), clusters)
    # All original indices present
    all_idx = sort(vcat(clusters...))
    @test all_idx == [1, 2, 3, 4]

    # Single point
    clusters_single = default_cluster([[0.0]], [1])
    @test length(clusters_single) == 1
    @test clusters_single[1] == [1]

    # Empty mapper
    mg_empty = mapper_graph(Vector{Vector{Float64}}(), p -> p[1], 3, 0.3, default_cluster)
    @test num_nodes(mg_empty) == 0

    # :pca filter
    pts_2d = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0], [0.0, 2.0]]
    mg_pca = mapper_graph(pts_2d, :pca, 2, 0.3, default_cluster)
    @test num_nodes(mg_pca) >= 1

    # Show method
    @test occursin("MapperGraph", sprint(show, mg))

    # Short-form mapper_graph (default overlap/cluster)
    mg_short = mapper_graph(points, p -> p[1], 3)
    @test num_nodes(mg_short) >= 1
end
