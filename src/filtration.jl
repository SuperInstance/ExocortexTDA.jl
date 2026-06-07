# ============================================================
# filtration.jl — Filtration and Vietoris-Rips
# ============================================================

"""
    Filtration

A filtration: an ordered sequence of simplices with associated filtration values.
Simplices are sorted by filtration value, then by dimension, then lexicographically.
"""
struct Filtration
    simplices::Vector{Simplex}
    values::Vector{Float64}

    function Filtration(simplices_in::Vector{Simplex}, values_in::Vector{Float64})
        @assert length(simplices_in) == length(values_in) "Simplices and values must have same length"
        perm = sortperm(collect(1:length(simplices_in)), by=i -> (values_in[i], dim(simplices_in[i]), vertices(simplices_in[i])))
        new(simplices_in[perm], values_in[perm])
    end
end

"""Simplices with filtration value ≤ threshold."""
function simplices_below(f::Filtration, value::Float64)
    result = Simplex[]
    for i in 1:length(f.simplices)
        if f.values[i] <= value
            push!(result, f.simplices[i])
        end
    end
    return result
end

"""Build SimplicialComplex at a given threshold value."""
function complex_at(f::Filtration, value::Float64)
    sc = SimplicialComplex()
    for s in simplices_below(f, value)
        if !(has_simplex(sc, s))
            push!(sc.simplices, s)
            d = dim(s)
            if !haskey(sc.by_dim, d)
                sc.by_dim[d] = Set{Simplex}()
            end
            push!(sc.by_dim[d], s)
        end
    end
    return sc
end

import Base: length
"""Number of simplices in filtration."""
length(f::Filtration) = length(f.simplices)

"""Filtration value of the i-th simplex."""
filtration_value(f::Filtration, i::Int) = f.values[i]

"""
    VietorisRipsFiltration(points, max_dim, max_radius)

Build a Vietoris-Rips filtration from a point cloud.

Iterative algorithm:
1. Compute all pairwise distances.
2. Add all vertices at value 0.0.
3. For each pair, add edge at distance.
4. For each k-simplex (k ≥ 2), add at the maximum pairwise distance among its vertices.
"""
function VietorisRipsFiltration(points::Vector{Vector{Float64}}, max_dim::Int, max_radius::Float64)
    n = length(points)
    n == 0 && return Filtration(Simplex[], Float64[])

    # Precompute pairwise distances
    dist = Matrix{Float64}(undef, n, n)
    for i in 1:n
        dist[i, i] = 0.0
        for j in (i + 1):n
            d = euclidean_dist(points[i], points[j])
            dist[i, j] = d
            dist[j, i] = d
        end
    end

    simplices_list = Simplex[]
    values_list = Float64[]

    # Add vertices
    for i in 1:n
        push!(simplices_list, Simplex(i))
        push!(values_list, 0.0)
    end

    # Add edges
    for i in 1:n
        for j in (i + 1):n
            d = dist[i, j]
            if d <= max_radius
                push!(simplices_list, Simplex(i, j))
                push!(values_list, d)
            end
        end
    end

    # Add higher simplices iteratively
    for k in 2:max_dim
        # Collect existing (k-1)-simplices
        prev_simps = filter(s -> dim(s) == k - 1, simplices_list)
        # Try adding each vertex to form a k-simplex
        seen = Set{Vector{Int}}()
        for s in prev_simps
            v = vertices(s)
            for i in 1:n
                if !(i in v)
                    new_verts = sort(vcat(v, [i]))
                    if !(new_verts in seen)
                        # Check all edges exist and compute max distance
                        valid = true
                        max_d = 0.0
                        for a in 1:length(new_verts)
                            for b in (a + 1):length(new_verts)
                                d = dist[new_verts[a], new_verts[b]]
                                if d > max_radius
                                    valid = false
                                    break
                                end
                                if d > max_d
                                    max_d = d
                                end
                            end
                            valid || break
                        end
                        if valid
                            push!(seen, new_verts)
                            push!(simplices_list, Simplex(new_verts))
                            push!(values_list, max_d)
                        end
                    end
                end
            end
        end
    end

    return Filtration(simplices_list, values_list)
end

# Multiple dispatch: accept Matrix too
function VietorisRipsFiltration(points::Matrix{Float64}, max_dim::Int, max_radius::Float64)
    pts = [points[i, :] for i in 1:size(points, 1)]
    return VietorisRipsFiltration(pts, max_dim, max_radius)
end

# Multiple dispatch: default max_radius
VietorisRipsFiltration(points::Vector{Vector{Float64}}, max_dim::Int) =
    VietorisRipsFiltration(points, max_dim, Inf)

VietorisRipsFiltration(points::Matrix{Float64}, max_dim::Int) =
    VietorisRipsFiltration(points, max_dim, Inf)

"""Euclidean distance between two vectors."""
function euclidean_dist(a::Vector{Float64}, b::Vector{Float64})
    s = 0.0
    for i in 1:length(a)
        d = a[i] - b[i]
        s += d * d
    end
    return sqrt(s)
end

euclidean_dist(a::Vector{<:Real}, b::Vector{<:Real}) = euclidean_dist(Float64.(a), Float64.(b))
