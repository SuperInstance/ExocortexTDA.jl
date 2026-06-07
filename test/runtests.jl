using ExocortexTDA
using Test

@testset "ExocortexTDA.jl" begin
    include("simplex_test.jl")
    include("filtration_test.jl")
    include("boundary_test.jl")
    include("barcode_test.jl")
    include("betti_test.jl")
    include("mapper_test.jl")
    include("distance_test.jl")
end
