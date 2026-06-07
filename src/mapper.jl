# ============================================================
# mapper.jl — Mapper algorithm
# ============================================================

"""
    MapperGraph

Result of the Mapper algorithm. Nodes are clusters of data points,
edges connect clusters with non-empty overlap.
"""
struct MapperGraph
    nodes::Vector{Set{Int}}       # Each node is a set of point indices
    edges::Set{Tuple{Int, Int}}   # Edges between nodes (sorted pairs)
    filter_values::Vector{Float64} # Filter value for each node (centroid)
end

"""
    mapper_graph(points, filter_func, n_intervals, overlap, cluster_func)

Mapper algorithm: filter → cover → cluster → nerve.

Multiple dispatch on `filter_func`:
- `Function` → apply directly to each point
- `Symbol` → use built-in (:pca, :eccentricity, :distance_to_mean)
"""
function mapper_graph(points::Vector{Vector{Float64}},
                      filter_func::Function,
                      n_intervals::Int,
                      overlap::Float64,
                      cluster_func::Function)
    n = length(points)
    n == 0 && return MapperGraph(Set{Int}[], Set{Tuple{Int,Int}}(), Float64[])

    # Step 1: Compute filter values
    fvals = [filter_func(p) for p in points]

    # Step 2: Create overlapping intervals
    f_min, f_max = extrema(fvals)
    f_range = f_max - f_min
    f_range == 0.0 && (f_range = 1.0)  # degenerate case

    step = f_range / n_intervals
    interval_len = step * (1.0 + overlap)

    intervals = Tuple{Float64, Float64}[]
    for i in 0:(n_intervals - 1)
        lo = f_min + i * step - overlap * step / 2
        hi = lo + interval_len
        push!(intervals, (lo, hi))
    end

    # Step 3: For each interval, cluster the points whose filter value falls in it
    all_nodes = Set{Int}[]
    all_fvals = Float64[]
    node_intervals = Vector{Int}[]  # which intervals each point belongs to

    for point_idx in 1:n
        push!(node_intervals, Int[])
    end

    for (int_idx, (lo, hi)) in enumerate(intervals)
        # Find points in this interval
        point_indices = Int[]
        for i in 1:n
            if fvals[i] >= lo && fvals[i] <= hi
                push!(point_indices, i)
            end
        end

        isempty(point_indices) && continue

        # Cluster these points
        clusters = cluster_func(points, point_indices)

        for cluster_indices in clusters
            node_id = length(all_nodes) + 1
            push!(all_nodes, Set(cluster_indices))
            # Average filter value
            avg_f = sum(fvals[i] for i in cluster_indices) / length(cluster_indices)
            push!(all_fvals, avg_f)
            for pi in cluster_indices
                push!(node_intervals[pi], node_id)
            end
        end
    end

    # Step 4: Build nerve — connect nodes with overlapping points
    edges = Set{Tuple{Int, Int}}()
    for i in 1:length(all_nodes)
        for j in (i + 1):length(all_nodes)
            if !isempty(all_nodes[i] ∩ all_nodes[j])
                push!(edges, (i, j))
            end
        end
    end

    return MapperGraph(all_nodes, edges, all_fvals)
end

# Multiple dispatch: Symbol filter functions
function mapper_graph(points::Vector{Vector{Float64}},
                      filter_sym::Symbol,
                      n_intervals::Int,
                      overlap::Float64,
                      cluster_func::Function)
    ff = _builtin_filter(filter_sym, points)
    return mapper_graph(points, ff, n_intervals, overlap, cluster_func)
end

# Multiple dispatch: Matrix points
function mapper_graph(points::Matrix{Float64},
                      filter_func,
                      n_intervals::Int,
                      overlap::Float64,
                      cluster_func::Function)
    pts = [points[i, :] for i in 1:size(points, 1)]
    return mapper_graph(pts, filter_func, n_intervals, overlap, cluster_func)
end

# Multiple dispatch: default overlap and cluster
mapper_graph(points::Vector{Vector{Float64}}, filter_func, n_intervals::Int) =
    mapper_graph(points, filter_func, n_intervals, 0.3, default_cluster)

mapper_graph(points::Matrix{Float64}, filter_func, n_intervals::Int) =
    mapper_graph(points, filter_func, n_intervals, 0.3, default_cluster)

"""
    _builtin_filter(sym, points) -> Function

Return a built-in filter function.
"""
function _builtin_filter(sym::Symbol, points::Vector{Vector{Float64}})
    if sym == :pca
        # First principal component (approximation via max-variance direction)
        centroid = _centroid(points)
        centered = [p .- centroid for p in points]
        # Power iteration for first PC
        d = length(points[1])
        v = randn(d)
        v /= sqrt(sum(v .^ 2))
        for _ in 1:50
            new_v = zeros(d)
            for p in centered
                proj = sum(p .* v)
                new_v .+= proj .* p
            end
            norm_v = sqrt(sum(new_v .^ 2))
            norm_v == 0.0 && break
            v = new_v / norm_v
        end
        return p -> sum((p .- centroid) .* v)
    elseif sym == :eccentricity
        # Eccentricity: average distance to all other points
        return p -> begin
            s = 0.0
            for q in points
                s += sqrt(sum((p .- q) .^ 2))
            end
            return s / length(points)
        end
    elseif sym == :distance_to_mean
        centroid = _centroid(points)
        return p -> sqrt(sum((p .- centroid) .^ 2))
    else
        error("Unknown filter: $sym. Use :pca, :eccentricity, or :distance_to_mean")
    end
end

function _centroid(points::Vector{Vector{Float64}})
    d = length(points[1])
    c = zeros(d)
    for p in points
        c .+= p
    end
    return c / length(points)
end

"""
    default_cluster(points, indices) -> Vector{Vector{Int}}

Simple single-linkage clustering with a distance threshold.
Splits points into clusters using a greedy connected-components approach.
"""
function default_cluster(points::Vector{Vector{Float64}}, indices::Vector{Int})
    isempty(indices) && return Vector{Int}[]
    length(indices) == 1 && return [indices]

    # Compute average pairwise distance for threshold
    total_dist = 0.0
    count = 0
    for i in 1:length(indices)
        for j in (i + 1):length(indices)
            total_dist += sqrt(sum((points[indices[i]] .- points[indices[j]]) .^ 2))
            count += 1
        end
    end
    threshold = count > 0 ? total_dist / count * 1.5 : 0.0

    # Union-Find
    parent = Dict(i => i for i in indices)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end

    for i in 1:length(indices)
        for j in (i + 1):length(indices)
            d = sqrt(sum((points[indices[i]] .- points[indices[j]]) .^ 2))
            if d <= threshold
                ri, rj = find(indices[i]), find(indices[j])
                if ri != rj
                    parent[ri] = rj
                end
            end
        end
    end

    # Collect clusters
    clusters = Dict{Int, Vector{Int}}()
    for idx in indices
        root = find(idx)
        if !haskey(clusters, root)
            clusters[root] = Int[]
        end
        push!(clusters[root], idx)
    end

    return collect(values(clusters))
end

"""Number of nodes in Mapper graph."""
num_nodes(mg::MapperGraph) = length(mg.nodes)

"""Number of edges in Mapper graph."""
num_edges(mg::MapperGraph) = length(mg.edges)

import Base: show
show(io::IO, mg::MapperGraph) = print(io, "MapperGraph($(num_nodes(mg)) nodes, $(num_edges(mg)) edges)")
