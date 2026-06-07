using ExocortexTDA
using Test

@testset "Simplex" begin
    # 0-simplex (vertex)
    v = Simplex(3)
    @test dim(v) == 0
    @test vertices(v) == [3]

    # 1-simplex (edge)
    e = Simplex(1, 3)
    @test dim(e) == 1
    @test vertices(e) == [1, 3]

    # Edge auto-sorts
    e2 = Simplex(5, 2)
    @test vertices(e2) == [2, 5]

    # General simplex
    t = Simplex([1, 2, 3])
    @test dim(t) == 2
    @test vertices(t) == [1, 2, 3]

    # Equality and hashing
    @test Simplex(1, 2) == Simplex(2, 1)
    @test hash(Simplex(1, 2)) == hash(Simplex(2, 1))

    # Boundary of a vertex is empty
    @test isempty(boundary(Simplex(1)))

    # Boundary of an edge: two vertices with signs
    bnd = boundary(Simplex(1, 2))
    @test length(bnd) == 2
    for p in bnd
        @test dim(p[1]) == 0
    end

    # Boundary of a triangle: three edges with signs
    bnd_t = boundary(Simplex([1, 2, 3]))
    @test length(bnd_t) == 3
    for p in bnd_t
        @test dim(p[1]) == 1
    end
    # Signs alternate: +1, -1, +1
    @test bnd_t[1][2] == 1
    @test bnd_t[2][2] == -1
    @test bnd_t[3][2] == 1

    # Mod-2 boundary
    @test length(boundary_mod2(Simplex([1, 2, 3]))) == 3

    # SimplicialComplex
    sc = SimplicialComplex()
    add_simplex!(sc, Simplex([1, 2, 3]))
    @test has_simplex(sc, Simplex(1, 2))
    @test has_simplex(sc, Simplex(1))
    @test length(sc) == 7  # 3 vertices + 3 edges + 1 triangle

    # f-vector
    @test f_vector(sc) == [3, 3, 1]

    # Euler characteristic of a single triangle: 3-3+1 = 1
    @test euler_characteristic(sc) == 1
end
