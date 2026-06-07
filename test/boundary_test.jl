using ExocortexTDA
using Test

@testset "Boundary" begin
    # Single edge: boundary matrix should have one column with one entry
    sc = SimplicialComplex()
    add_simplex!(sc, Simplex(1, 2))
    M = boundary_matrix(sc)
    @test length(M.columns) == 3  # 2 vertices + 1 edge

    # The edge column should have entries pointing to the two vertex rows
    edge_idx = findfirst(s -> dim(s) == 1, M.simplex_order)
    @test edge_idx !== nothing
    @test length(M.columns[edge_idx]) == 2

    # Vertex columns should be empty
    for i in 1:length(M.columns)
        if dim(M.simplex_order[i]) == 0
            @test isempty(M.columns[i])
        end
    end

    # low function
    @test low(M, edge_idx) > 0
    for i in 1:length(M.columns)
        if dim(M.simplex_order[i]) == 0
            @test low(M, i) == 0
        end
    end

    # Reduction of a single triangle
    sc2 = SimplicialComplex()
    add_simplex!(sc2, Simplex([1, 2, 3]))
    M2 = boundary_matrix(sc2)
    pairs = reduce_matrix!(M2)

    # Triangle has 7 simplices: 3 vertices, 3 edges, 1 triangle
    @test length(M2.columns) == 7
    # After reduction, should get persistence pairs
    @test length(pairs) >= 1

    # Empty complex
    sc_empty = SimplicialComplex()
    M_empty = boundary_matrix(sc_empty)
    @test length(M_empty.columns) == 0

    # Two disconnected edges
    sc3 = SimplicialComplex()
    add_simplex!(sc3, Simplex(1, 2))
    add_simplex!(sc3, Simplex(3, 4))
    M3 = boundary_matrix(sc3)
    pairs3 = reduce_matrix!(M3)
    @test length(M3.columns) == 6  # 4 vertices + 2 edges

    # After reduction: each edge pairs with one of its vertices,
    # leaving 2 unpaired vertices (H_0 generators)
    unpaired = unpaired_columns(M3, pairs3)
    @test length(unpaired) == 2
end
