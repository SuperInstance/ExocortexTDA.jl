# ============================================================
# distance.jl — Bottleneck and Wasserstein distances
# ============================================================

"""
    bottleneck_distance(d1::PersistenceDiagram, d2::PersistenceDiagram) -> Float64

Compute the bottleneck distance between two persistence diagrams.
Uses greedy matching for efficiency.
"""
function bottleneck_distance(d1::PersistenceDiagram, d2::PersistenceDiagram)
    isempty(d1.points) && isempty(d2.points) && return 0.0
    isempty(d1.points) && return _max_diag_distance(d2)
    isempty(d2.points) && return _max_diag_distance(d1)

    pts1 = collect(d1.points)
    pts2 = collect(d2.points)
    n1, n2 = length(pts1), length(pts2)

    # Build all possible assignments sorted by distance
    assignments = Tuple{Float64, Int, Int}[]
    for i in 1:n1
        for j in 1:n2
            d = linf_dist(pts1[i], pts2[j])
            push!(assignments, (d, i, j))
        end
        # Distance to diagonal
        d = diag_dist(pts1[i])
        push!(assignments, (d, i, -1))
    end
    for j in 1:n2
        d = diag_dist(pts2[j])
        push!(assignments, (d, -1, j))
    end

    sort!(assignments, by=x -> x[1])

    used1 = falses(n1)
    used2 = falses(n2)
    max_dist = 0.0

    for (d, i, j) in assignments
        if i > 0 && j > 0
            if !used1[i] && !used2[j]
                used1[i] = true
                used2[j] = true
                max_dist = max(max_dist, d)
            end
        elseif i > 0
            if !used1[i]
                used1[i] = true
                max_dist = max(max_dist, d)
            end
        elseif j > 0
            if !used2[j]
                used2[j] = true
                max_dist = max(max_dist, d)
            end
        end
    end

    return max_dist
end

"""
    wasserstein_distance(d1::PersistenceDiagram, d2::PersistenceDiagram, p::Int) -> Float64

Compute W_p Wasserstein distance using greedy matching.
"""
function wasserstein_distance(d1::PersistenceDiagram, d2::PersistenceDiagram, p::Int)
    isempty(d1.points) && isempty(d2.points) && return 0.0
    isempty(d1.points) && return _wasserstein_diag(d2, p)
    isempty(d2.points) && return _wasserstein_diag(d1, p)

    pts1 = collect(d1.points)
    pts2 = collect(d2.points)
    n1, n2 = length(pts1), length(pts2)

    # Greedy matching
    used1 = falses(n1)
    used2 = falses(n2)
    total = 0.0

    # Build sorted assignment list
    assignments = Tuple{Float64, Int, Int}[]
    for i in 1:n1
        for j in 1:n2
            d = linf_dist(pts1[i], pts2[j])
            push!(assignments, (d, i, j))
        end
    end
    sort!(assignments, by=x -> x[1])

    # Greedy assign
    for (d, i, j) in assignments
        if !used1[i] && !used2[j]
            total += d^p
            used1[i] = true
            used2[j] = true
        end
    end

    # Unmatched points → diagonal distance
    for i in 1:n1
        if !used1[i]
            total += diag_dist(pts1[i])^p
        end
    end
    for j in 1:n2
        if !used2[j]
            total += diag_dist(pts2[j])^p
        end
    end

    return total^(1.0 / p)
end

# Default p=2
wasserstein_distance(d1::PersistenceDiagram, d2::PersistenceDiagram) =
    wasserstein_distance(d1, d2, 2)

"""L∞ distance between two persistence points."""
function linf_dist(p1::Tuple{Float64,Float64}, p2::Tuple{Float64,Float64})
    b1, d1 = p1
    b2, d2 = p2
    db = abs(b1 - b2)
    if isinf(d1) && isinf(d2)
        return db
    elseif isinf(d1) || isinf(d2)
        return Inf
    else
        dd = abs(d1 - d2)
        return max(db, dd)
    end
end

"""Distance from a persistence point to the diagonal."""
function diag_dist(p::Tuple{Float64,Float64})
    b, d = p
    if isinf(d)
        return 0.0  # Essential cycles live on the diagonal at infinity
    end
    return abs(d - b) / 2.0
end

"""Max diagonal distance for all points in a diagram."""
function _max_diag_distance(d::PersistenceDiagram)
    isempty(d.points) && return 0.0
    return maximum(diag_dist(p) for p in d.points)
end

function _wasserstein_diag(d::PersistenceDiagram, p::Int)
    total = 0.0
    for pt in d.points
        total += diag_dist(pt)^p
    end
    return total^(1.0 / p)
end
