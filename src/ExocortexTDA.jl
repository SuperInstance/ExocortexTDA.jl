# ============================================================
# ExocortexTDA.jl — Topological Data Analysis for the Exocortex
# ============================================================

module ExocortexTDA

# Simplex and simplicial complex
export Simplex, SimplicialComplex,
       dim, vertices, boundary, boundary_mod2,
       add_simplex!, simplices, faces, has_simplex,
       max_dim, f_vector, euler_characteristic

# Filtration
export Filtration, VietorisRipsFiltration,
       complex_at, simplices_below, filtration_value, euclidean_dist

# Boundary matrix
export BoundaryMatrix, boundary_matrix, reduce_matrix!, low, unpaired_columns

# Barcode
export PersistenceBarcode, PersistenceDiagram,
       num_bars, num_points, bars_in_dim,
       compute_persistence, plot_ascii

# Betti
export betti_numbers, betti_curve, print_betti

# Mapper
export MapperGraph, mapper_graph, num_nodes, num_edges, default_cluster

# Distance
export bottleneck_distance, wasserstein_distance, linf_dist, diag_dist

include("simplex.jl")
include("filtration.jl")
include("boundary.jl")
include("barcode.jl")
include("betti.jl")
include("mapper.jl")
include("distance.jl")

end # module ExocortexTDA
