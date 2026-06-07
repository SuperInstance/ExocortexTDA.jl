# ExocortexTDA.jl

Topological data analysis for the exocortex — proving multiple dispatch makes TDA readable and composable.

[![Julia](https://img.shields.io/badge/Julia-1.10+-9558B2)](https://julialang.org)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-zero-green)]()

```
ExocortexTDA.jl
├── src/
│   ├── simplex.jl      Simplex + SimplicialComplex (multiple dispatch constructors)
│   ├── filtration.jl    Filtration + Vietoris-Rips (multiple dispatch input types)
│   ├── boundary.jl      Boundary matrix reduction (mod-2, iterative)
│   ├── barcode.jl       PersistenceBarcode + PersistenceDiagram + ASCII plots
│   ├── betti.jl         Betti numbers + Betti curves (4 dispatch variants)
│   ├── mapper.jl        Mapper algorithm (Symbol/Function/Matrix dispatch)
│   └── distance.jl      Bottleneck + Wasserstein distances
├── test/                50+ tests covering every module
├── Project.toml
└── README.md
```

## Why Julia for TDA?

> *Multiple dispatch makes the topology match the mathematics.*

In TDA, the same operation has different meanings depending on context. Computing Betti numbers from a barcode, a simplicial complex, or a filtration are conceptually identical but algorithmically distinct. Julia's multiple dispatch lets you write `betti_numbers(x)` and get the right implementation automatically — no visitor patterns, no type tags, no switch statements.

```julia
betti_numbers(barcode::PersistenceBarcode, threshold)  # count alive bars
betti_numbers(sc::SimplicialComplex)                    # reduce boundary matrix
betti_numbers(f::Filtration, value)                     # compute persistence first
```

The same principle applies to simplex constructors (`Simplex(3)` → vertex, `Simplex(1,2)` → edge, `Simplex([1,2,3])` → triangle), filtration inputs (`Vector{Vector}` or `Matrix`), mapper filters (`Function` or `:pca`), and persistence computation (`Filtration` or `SimplicialComplex`).

## Installation

```julia
# In Julia REPL
using Pkg
Pkg.add(url="https://github.com/SuperInstance/ExocortexTDA.jl")
```

Or clone and develop:

```bash
git clone https://github.com/SuperInstance/ExocortexTDA.jl
cd ExocortexTDA.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Quick Start

```julia
using ExocortexTDA

# Build a Vietoris-Rips filtration from a point cloud
points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
f = VietorisRipsFiltration(points, 2, 3.0)

# Compute persistent homology
bc = compute_persistence(f)

# Visualize the barcode
println(plot_ascii(bc))

# Betti numbers at threshold 1.5
bn = betti_numbers(bc, 1.5)
println("β₀ = $(get(bn, 0, 0)), β₁ = $(get(bn, 1, 0))")

# Mapper graph
mg = mapper_graph(points, :pca, 3)
println(mg)
```

---

## Architecture

```
  ┌─────────────┐
  │  Point Cloud │  ── or any metric space
  └──────┬──────┘
         │  VietorisRipsFiltration()
         ▼
  ┌─────────────┐
  │  Filtration  │  sorted simplices + values
  └──────┬──────┘
         │  boundary_matrix()
         ▼
  ┌─────────────┐
  │   Boundary   │  mod-2 sparse matrix
  │    Matrix    │
  └──────┬──────┘
         │  reduce_matrix!()
         ▼
  ┌─────────────┐
  │  Persistence │  (birth, death) pairs
  │   Barcode    │
  └──┬───┬───┬──┘
     │   │   │
     ▼   ▼   ▼
  ┌────┐┌────┐┌─────────┐
  │Betti││Mapper││Distance │
  │ Nums││Graph││Bottleneck│
  └─────┘└─────┘│Wasserstein│
                └──────────┘
```

**Data flow:** raw points → filtration → boundary matrix → reduction → barcode → downstream analysis.

Each arrow is a function call with multiple dispatch variants. The pipeline is composable: you can enter at any stage with the appropriate data structure.

---

## Mathematical Background

### Topological Data Analysis

Topological data analysis (TDA) applies algebraic topology to study the "shape" of data. The central insight is that data has topology — holes, voids, connected components — that persist across scales and capture structural features invisible to traditional statistics.

The key object is a **simplicial complex** K, a collection of simplices (points, edges, triangles, tetrahedra, ...) closed under taking faces. Simplicial complexes generalize graphs to higher dimensions.

**Homology** H_k(K) captures k-dimensional topological features:
- H_0: connected components
- H_1: loops (1-dimensional holes)
- H_2: voids (2-dimensional cavities)
- H_k: k-dimensional "holes"

Formally, H_k(X) = ker(∂_k) / im(∂_{k+1}), where ∂_k is the boundary operator mapping k-chains to (k-1)-chains. A cycle in ker(∂_k) that is not a boundary (not in im(∂_{k+1})) represents a genuine hole.

### Persistent Homology

A single simplicial complex gives a snapshot of topology at one scale. **Persistent homology** tracks how topology evolves across a filtered sequence of complexes:

K_0 ⊆ K_1 ⊆ K_2 ⊆ ... ⊆ K_n

The output is a **persistence barcode**: a multiset of intervals [b_i, d_i) where b_i is the filtration value where a feature is born and d_i is where it dies. Long bars represent significant topological features; short bars are noise.

Equivalently, the barcode can be represented as a **persistence diagram**: a multiset of points (b_i, d_i) in the plane above the diagonal. Points far from the diagonal correspond to persistent features.

The boundary matrix reduction algorithm computes persistent homology by:
1. Constructing the mod-2 boundary matrix from the filtered complex
2. Column-reducing it (iterative left-to-right Gaussian elimination over ℤ₂)
3. Reading off persistence pairs from the reduced matrix

### The Mapper Algorithm

The Mapper algorithm (Singh, Mémoli & Carlsson, 2007) constructs a combinatorial graph summary of high-dimensional data:

1. **Filter:** Apply a function f: X → ℝ to the data
2. **Cover:** Partition the range of f into overlapping intervals
3. **Cluster:** Within each interval, cluster the pre-image
4. **Nerve:** Connect clusters with non-empty overlap

The resulting Mapper graph captures the topology of the data in a combinatorial form amenable to visualization and analysis.

### Distance Metrics

Two standard metrics compare persistence diagrams:

- **Bottleneck distance:** d_∞(D₁, D₂) = min_η max_x ‖x - η(x)‖_∞, where η ranges over bijections from D₁ ∪ Δ to D₂ ∪ Δ (Δ = diagonal). The bottleneck distance captures the single worst-matched feature.

- **Wasserstein distance:** W_p(D₁, D₂) = (min_η Σ ‖x - η(x)‖_∞^p)^{1/p}. The Wasserstein distance aggregates all matching costs, giving a more holistic comparison.

Both satisfy the stability theorem: small perturbations of the input produce small changes in the diagram distance.

---

## API Reference

### Simplex and SimplicialComplex

```julia
Simplex(3)              # 0-simplex (vertex)
Simplex(1, 2)           # 1-simplex (edge), auto-sorted
Simplex([1, 2, 3])      # 2-simplex (triangle)

dim(s)                  # dimension (0, 1, 2, ...)
vertices(s)             # sorted vertex list
boundary(s)             # [(face, sign), ...]
boundary_mod2(s)        # [face1, face2, ...]

sc = SimplicialComplex()
add_simplex!(sc, Simplex([1, 2, 3]))
simplices(sc)           # all simplices
simplices(sc, 1)        # edges only
has_simplex(sc, Simplex(1, 2))
f_vector(sc)            # [n₀, n₁, n₂, ...]
euler_characteristic(sc) # Σ(-1)^k f_k
```

### Filtration

```julia
# From explicit simplex-value pairs
f = Filtration([Simplex(1), Simplex(1,2)], [0.0, 1.5])

# Vietoris-Rips from point cloud
f = VietorisRipsFiltration(points, max_dim, max_radius)
f = VietorisRipsFiltration(points, max_dim)  # unlimited radius

complex_at(f, 0.5)      # SimplicialComplex at threshold
simplices_below(f, 0.5) # simplices with value ≤ threshold
```

### Boundary Matrix

```julia
M = boundary_matrix(sc)  # from SimplicialComplex
M = boundary_matrix(f)    # from Filtration

pairs = reduce_matrix!(M)  # [(birth_col, death_col), ...]
low(M, j)                   # lowest 1 in column j
unpaired_columns(M, pairs)  # essential cycles
```

### Barcode and Diagram

```julia
bc = compute_persistence(f)  # from Filtration
bc = compute_persistence(sc)  # from SimplicialComplex

num_bars(bc)             # total
num_bars(bc, 1)          # in dimension 1
bars_in_dim(bc, 0)       # intervals for H₀

pd = PersistenceDiagram(bc)
println(plot_ascii(bc))  # ASCII barcode
println(plot_ascii(pd))  # ASCII diagram
```

### Betti Numbers

```julia
betti_numbers(barcode, threshold)    # Dict(dim => count)
betti_numbers(sc)                    # from SimplicialComplex
betti_numbers(f, value)              # from Filtration

betti_curve(barcode, thresholds)     # Dict(dim => Vector{Int})
betti_curve(f, thresholds)           # from Filtration
```

### Mapper

```julia
# Function filter
mg = mapper_graph(points, p -> p[1], n_intervals, overlap, cluster_func)

# Built-in filter (Symbol)
mg = mapper_graph(points, :pca, n_intervals, overlap, cluster_func)
# :pca, :eccentricity, :distance_to_mean

# Matrix input
mg = mapper_graph(matrix, filter, n_intervals)

# Default overlap and clustering
mg = mapper_graph(points, :pca, 5)

num_nodes(mg)  # node count
num_edges(mg)  # edge count
```

### Distance Metrics

```julia
bottleneck_distance(d1, d2)           # d_∞
wasserstein_distance(d1, d2)          # W₂ (default)
wasserstein_distance(d1, d2, 1)       # W₁
```

---

## Examples

### Example 1: Persistent Homology of a Point Cloud

```julia
using ExocortexTDA

# Generate points on a circle (should have H₀=1, H₁=1)
n_points = 20
angles = [2π * i / n_points for i in 0:(n_points-1)]
circle = [[cos(a), sin(a)] for a in angles]

# Build Vietoris-Rips filtration
f = VietorisRipsFiltration(circle, 2, 2.0)
println("Filtration has $(length(f)) simplices")

# Compute persistent homology
bc = compute_persistence(f)
println("\nPersistence barcode:")
println(plot_ascii(bc))

# Betti numbers at various scales
for t in [0.1, 0.5, 1.0, 1.5]
    bn = betti_numbers(bc, t)
    println("t=$t: β₀=$(get(bn, 0, 0)), β₁=$(get(bn, 1, 0))")
end

# Betti curve
thresholds = collect(0.0:0.1:2.0)
curve = betti_curve(bc, thresholds)
println("\nH₀ curve: $(curve[0])")
if haskey(curve, 1)
    println("H₁ curve: $(curve[1])")
end
```

### Example 2: Mapper Graph with PCA Filter

```julia
using ExocortexTDA

# Generate two Gaussian clusters
cluster1 = [[randn() * 0.1 - 2.0, randn() * 0.1] for _ in 1:30]
cluster2 = [[randn() * 0.1 + 2.0, randn() * 0.1] for _ in 1:30]
points = vcat(cluster1, cluster2)

# Mapper with PCA filter
mg = mapper_graph(points, :pca, 4, 0.3, default_cluster)
println("Mapper graph: $(num_nodes(mg)) nodes, $(num_edges(mg)) edges")

# Print each node
for (i, node) in enumerate(mg.nodes)
    println("  Node $i: $(length(node)) points, filter=$(round(mg.filter_values[i], digits=2))")
end

# Print edges
println("\nEdges:")
for (i, j) in mg.edges
    overlap = length(mg.nodes[i] ∩ mg.nodes[j])
    println("  $i ↔ $j (overlap: $overlap points)")
end

# Also try eccentricity filter
mg2 = mapper_graph(points, :eccentricity, 3, 0.3, default_cluster)
println("\nWith eccentricity: $(num_nodes(mg2)) nodes, $(num_edges(mg2)) edges")
```

### Example 3: Betti Curve Analysis

```julia
using ExocortexTDA

# Point cloud sampling from a torus (2D projection)
n = 50
torus_pts = [[(2 + cos(2π*i/n))*cos(2π*j/n),
              (2 + cos(2π*i/n))*sin(2π*j/n)]
             for i in 1:n, j in 1:5]

points = vec(torus_pts)

# Build filtration
f = VietorisRipsFiltration(points, 2, 1.5)

# Compute persistence
bc = compute_persistence(f)

# ASCII barcode
println("Barcode:")
println(plot_ascii(bc))

# Betti curve at fine resolution
thresholds = collect(0.0:0.05:1.5)
curve = betti_curve(bc, thresholds)

println("\nBetti curve (H₀):")
for (i, t) in enumerate(thresholds)
    bars = get(curve, 0, zeros(Int, length(thresholds)))
    if i % 5 == 1
        bar_str = "█" ^ bars[i]
        println("  t=$(round(t, digits=2)): $bar_str ($(bars[i]))")
    end
end

# Direct Betti numbers from filtration
for t in [0.2, 0.5, 1.0]
    bn = betti_numbers(f, t)
    println("\nAt t=$t:")
    for d in sort(collect(keys(bn)))
        println("  β_$d = $(bn[d])")
    end
end
```

### Example 4: Comparing Two Datasets via Bottleneck Distance

```julia
using ExocortexTDA

# Dataset 1: circle
n = 15
circle = [[cos(2π*i/n), sin(2π*i/n)] for i in 0:(n-1)]
f1 = VietorisRipsFiltration(circle, 2, 2.0)
bc1 = compute_persistence(f1)
pd1 = PersistenceDiagram(bc1)

# Dataset 2: circle with noise
noisy_circle = [[cos(2π*i/n) + 0.1*randn(), sin(2π*i/n) + 0.1*randn()] for i in 0:(n-1)]
f2 = VietorisRipsFiltration(noisy_circle, 2, 2.0)
bc2 = compute_persistence(f2)
pd2 = PersistenceDiagram(bc2)

# Dataset 3: line (different topology)
line_pts = [[i/n, 0.0] for i in 0:n]
f3 = VietorisRipsFiltration(line_pts, 2, 2.0)
bc3 = compute_persistence(f3)
pd3 = PersistenceDiagram(bc3)

# Compare
bn_circle_noisy = bottleneck_distance(pd1, pd2)
bn_circle_line = bottleneck_distance(pd1, pd3)
bn_noisy_line = bottleneck_distance(pd2, pd3)

println("Bottleneck distances:")
println("  Circle ↔ Noisy circle: $(round(bn_circle_noisy, digits=4))")
println("  Circle ↔ Line:         $(round(bn_circle_line, digits=4))")
println("  Noisy circle ↔ Line:   $(round(bn_noisy_line, digits=4))")

# Wasserstein distances (p=1 and p=2)
println("\nWasserstein distances:")
w1_cn = wasserstein_distance(pd1, pd2, 1)
w2_cn = wasserstein_distance(pd1, pd2, 2)
println("  Circle ↔ Noisy: W₁=$(round(w1_cn, digits=4)), W₂=$(round(w2_cn, digits=4))")

w1_cl = wasserstein_distance(pd1, pd3, 1)
w2_cl = wasserstein_distance(pd1, pd3, 2)
println("  Circle ↔ Line:  W₁=$(round(w1_cl, digits=4)), W₂=$(round(w2_cl, digits=4))")
```

---

## Design Decisions

### Why Multiple Dispatch?

TDA naturally involves operations that vary by input type. `betti_numbers` means:
- Count alive bars in a barcode (trivial)
- Reduce a boundary matrix from a complex (moderate)
- Compute persistence then count (expensive)

Multiple dispatch selects the right algorithm without runtime checks, producing type-stable, JIT-friendly code. The compiler sees each method as a separate function, enabling aggressive optimization.

### Why Mod-2 Coefficients?

We use ℤ₂ (mod-2) boundary matrices for simplicity and performance:
- No sign tracking needed (XOR replaces subtraction)
- Sufficient for most data analysis applications
- Matches the standard in computational topology (Edelsbrunner & Harer, 2010)

Extending to ℤ_p or ℤ coefficients would follow the same dispatch pattern.

### Why Iterative Algorithms?

All algorithms are iterative (no recursion):
- No stack overflow on large inputs
- Predictable memory usage
- Better cache behavior
- Easier to reason about complexity

The boundary reduction is the standard iterative left-to-right column algorithm.

### Why Greedy Matching?

For bottleneck and Wasserstein distances, we use greedy matching:
- O(n² log n) vs O(n³) for exact Hungarian algorithm
- Good approximation in practice
- Exact for bottleneck distance (by the augmenting path property)
- Zero external dependencies

For production use with large diagrams, one would want the Hopcroft-Karp algorithm or auction algorithm.

### Why Zero Dependencies?

ExocortexTDA.jl uses only Julia Base:
- No installation friction
- No version conflicts
- Compileable to static binaries
- Demonstrates Julia's standard library sufficiency

The tradeoff: performance-sensitive applications should use Ripser.jl or Eirene.jl for large-scale computations.

---

## Performance

Julia's JIT compilation makes ExocortexTDA.jl competitive for small-to-medium datasets:

- **Loop fusion:** Julia's `for` loops compile to tight machine code with no overhead
- **Cache-friendly arrays:** Simplices stored as sorted integer vectors for locality
- **Sparse representation:** Boundary matrix uses Set{Int} columns (no dense allocation)
- **Type stability:** Multiple dispatch ensures type-stable code paths
- **In-place reduction:** `reduce_matrix!` modifies the matrix in-place

For large point clouds (>10K points), use specialized tools:
- **Ripser.jl:** Optimized C++ bindings for Vietoris-Rips persistence
- **Eirene.jl:** Julia-native, handles larger complexes
- **PersistenceDiagrams.jl:** Efficient diagram operations

ExocortexTDA.jl prioritizes **readability and composability** over raw speed. The entire pipeline is Julia code you can read, modify, and extend.

---

## Comparison with Other TDA Tools

| Feature | ExocortexTDA.jl | Ripser.jl | Eirene.jl | scikit-tda (Python) | PersistenceDiagrams.jl |
|---|---|---|---|---|---|
| Language | Julia | Julia (C++ backend) | Julia | Python | Julia |
| Dependencies | 0 | Binary | Some | Many | Some |
| VR Persistence | ✓ (educational) | ✓ (optimized) | ✓ | ✓ | ✗ |
| Mapper | ✓ | ✗ | ✗ | ✓ | ✗ |
| Barcode plots | ASCII | ✗ | GUI | Matplotlib | Plots.jl |
| Multiple dispatch | Core design | Limited | Some | N/A | Some |
| Iterative algos | All | All | Mixed | Mixed | Mixed |
| Bottleneck dist | ✓ (greedy) | ✗ | ✓ | ✓ | ✓ |
| Readability | High | Low | Medium | Medium | Medium |

**ExocortexTDA.jl** is the teaching library. Use it to understand TDA, prototype ideas, and analyze small datasets. For production work, combine it with Ripser.jl for fast persistence computation and use ExocortexTDA for the downstream analysis (Betti curves, Mapper, distances).

---

## References

1. Edelsbrunner, H. & Harer, J. (2010). *Computational Topology: An Introduction*. American Mathematical Society. — The definitive textbook on persistent homology and the boundary matrix reduction algorithm.

2. Carlsson, G. (2009). "Topology and data." *Bulletin of the AMS*, 46(2), 255–308. — Foundational paper establishing TDA as a field, introducing the statistical perspective on persistent homology.

3. Singh, G., Mémoli, F. & Carlsson, G. (2007). "Topological methods for the analysis of high dimensional data sets and 3D object recognition." *Eurographics Symposium on Point-Based Graphics*. — Original Mapper algorithm paper.

4. Zomorodian, A. & Carlsson, G. (2005). "Computing persistent homology." *Discrete & Computational Geometry*, 33(2), 249–274. — Efficient algorithms for persistent homology computation.

5. Chazal, F. et al. (2013). "Persistence-based clustering in Riemannian manifolds." *Computational Geometry*, 46(3), 381–403. — Stability theorems for persistence-based methods.

6. Ghrist, R. (2008). "Barcodes: The persistent topology of data." *Bulletin of the AMS*, 45(1), 61–75. — Accessible survey of persistent homology with the barcode perspective.

7. Kerber, M., Morozov, D. & Nigmetov, A. (2017). "Geometry helps to compare persistence diagrams." *Journal of Experimental Algorithmics*, 22, 1–20. — Efficient algorithms for bottleneck and Wasserstein distances.

8. Bauer, U. (2021). "Ripser: efficient computation of Vietoris-Rips persistence barcodes." *Journal of Applied and Computational Topology*, 5, 391–423. — The Ripser algorithm and implementation that inspired our reduction approach.

---

## Glossary

| Term | Definition |
|---|---|
| **Simplex** | The convex hull of k+1 affinely independent points. A 0-simplex is a vertex, 1-simplex an edge, 2-simplex a triangle. |
| **Simplicial Complex** | A collection of simplices closed under taking faces. Generalizes graphs to higher dimensions. |
| **Boundary Operator** | ∂_k maps a k-simplex to the sum (mod 2) of its (k-1)-dimensional faces. Satisfies ∂² = 0. |
| **Homology** | H_k(X) = ker(∂_k) / im(∂_{k+1}). Measures k-dimensional topological features (holes). |
| **Betti Number** | β_k = rank(H_k). The number of independent k-dimensional holes. |
| **Filtration** | A nested sequence of simplicial complexes K_0 ⊆ K_1 ⊆ ... ⊆ K_n indexed by a parameter. |
| **Persistent Homology** | The study of how homology evolves across a filtration. Tracks birth and death of topological features. |
| **Persistence Barcode** | A multiset of intervals [b_i, d_i) representing the lifetime of each topological feature. |
| **Persistence Diagram** | A multiset of points (b_i, d_i) in the plane, equivalent to the barcode. |
| **Vietoris-Rips Complex** | A simplicial complex built from a point cloud: include a k-simplex when all pairwise distances ≤ ε. |
| **Mapper Graph** | A combinatorial graph summary of data: filter → cover → cluster → nerve. |
| **Bottleneck Distance** | d_∞(D₁, D₂) = min matching cost under L∞ norm. Captures the worst-matched feature. |
| **Wasserstein Distance** | W_p(D₁, D₂) = (Σ matching costs^p)^{1/p}. Aggregates all matching costs. |
| **Mod-2 Coefficients** | Arithmetic over ℤ₂ = {0, 1} with addition as XOR. Simplifies boundary computation. |
| **Nerve** | A simplicial complex built from a cover: k-simplices correspond to (k+1)-fold intersections of cover elements. |
| **F-vector** | [f₀, f₁, f₂, ...] where f_k = number of k-simplices in the complex. |
| **Euler Characteristic** | χ(K) = Σ(-1)^k f_k. A topological invariant related to Betti numbers by χ = Σ(-1)^k β_k. |
| **Column Reduction** | Left-to-right Gaussian elimination of the boundary matrix over ℤ₂. Produces the persistence pairing. |
| **Essential Cycle** | A homology class that persists to infinity (born but never dies). Corresponds to unpaired columns. |
| **Multiple Dispatch** | Julia's method of selecting function implementations based on the types of all arguments. |

---

## License

MIT

---

*Part of the [Exocortex](https://github.com/SuperInstance) project — topological tools for augmented cognition.*
