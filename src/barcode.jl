# ============================================================
# barcode.jl — Persistence barcodes and diagrams
# ============================================================

"""
    PersistenceBarcode

A collection of persistence intervals [birth, death).
Finite bars have death > birth; infinite bars have death = Inf.
"""
struct PersistenceBarcode
    intervals::Vector{Tuple{Float64, Float64}}
    dimensions::Vector{Int}
end

PersistenceBarcode(intervals::Vector{Tuple{Float64, Float64}}) =
    PersistenceBarcode(intervals, zeros(Int, length(intervals)))

"""Number of bars."""
num_bars(b::PersistenceBarcode) = length(b.intervals)

"""Number of bars in a given dimension."""
function num_bars(b::PersistenceBarcode, dim::Int)
    count = 0
    for i in 1:length(b.dimensions)
        if b.dimensions[i] == dim
            count += 1
        end
    end
    return count
end

"""Bars in a specific dimension."""
function bars_in_dim(b::PersistenceBarcode, d::Int)
    result = Tuple{Float64, Float64}[]
    for i in 1:length(b.intervals)
        if b.dimensions[i] == d
            push!(result, b.intervals[i])
        end
    end
    return result
end

"""
    PersistenceDiagram

Persistence diagram: multiset of points (birth, death) in the plane,
plus points at (birth, Inf) for essential cycles.
"""
struct PersistenceDiagram
    points::Vector{Tuple{Float64, Float64}}
    dimensions::Vector{Int}
end

PersistenceDiagram(points::Vector{Tuple{Float64, Float64}}) =
    PersistenceDiagram(points, zeros(Int, length(points)))

"""Convert barcode to persistence diagram."""
function PersistenceDiagram(b::PersistenceBarcode)
    PersistenceDiagram(copy(b.intervals), copy(b.dimensions))
end

"""Number of points in diagram."""
num_points(d::PersistenceDiagram) = length(d.points)

"""
    compute_persistence(f::Filtration) -> PersistenceBarcode

Compute persistent homology from a filtration.
Uses boundary matrix reduction internally.
"""
function compute_persistence(f::Filtration)
    M = boundary_matrix(f)
    pairs = reduce_matrix!(M)

    intervals = Tuple{Float64, Float64}[]
    dims = Int[]

    for (birth_idx, death_idx) in pairs
        birth_val = f.values[birth_idx]
        death_val = f.values[death_idx]
        # Dimension of the homology class
        d = dim(f.simplices[birth_idx])
        push!(intervals, (birth_val, death_val))
        push!(dims, d)
    end

    # Essential cycles (unpaired columns → infinite bars)
    unpaired = unpaired_columns(M, pairs)
    for idx in unpaired
        d = dim(f.simplices[idx])
        birth_val = f.values[idx]
        push!(intervals, (birth_val, Inf))
        push!(dims, d)
    end

    return PersistenceBarcode(intervals, dims)
end

# Multiple dispatch: from SimplicialComplex (single complex, trivial barcode)
function compute_persistence(sc::SimplicialComplex)
    M = boundary_matrix(sc)
    pairs = reduce_matrix!(M)

    intervals = Tuple{Float64, Float64}[]
    dims = Int[]
    all_simps = M.simplex_order

    for (b_idx, d_idx) in pairs
        d = dim(all_simps[b_idx])
        push!(intervals, (0.0, 1.0))
        push!(dims, d)
    end

    unpaired = unpaired_columns(M, pairs)
    for idx in unpaired
        d = dim(all_simps[idx])
        push!(intervals, (0.0, Inf))
        push!(dims, d)
    end

    return PersistenceBarcode(intervals, dims)
end

"""
    plot_ascii(barcode::PersistenceBarcode) -> String

ASCII visualization of a persistence barcode.
"""
function plot_ascii(barcode::PersistenceBarcode)
    isempty(barcode.intervals) && return "(empty barcode)"

    max_death = 0.0
    for (b, d) in barcode.intervals
        if isinf(d)
            max_death = max(max_death, b + 1.0)
        else
            max_death = max(max_death, d)
        end
    end

    width = 60
    lines = String[]

    for dim_val in sort(unique(barcode.dimensions))
        push!(lines, "─── Dimension $dim_val ───")
        dim_bars = filter(i -> barcode.dimensions[i] == dim_val, 1:length(barcode.intervals))
        for i in dim_bars
            b, d = barcode.intervals[i]
            d_str = isinf(d) ? "∞" : sprint(show, round(d, digits=2))
            b_pos = round(Int, b / max_death * (width - 1)) + 1
            b_pos = clamp(b_pos, 1, width)
            line = repeat(' ', b_pos - 1) * "|"
            if isinf(d)
                line *= repeat('─', width - b_pos) * "→ $b → ∞"
            else
                d_pos = round(Int, d / max_death * (width - 1)) + 1
                d_pos = clamp(d_pos, b_pos + 1, width)
                bar_len = d_pos - b_pos - 1
                line *= repeat('━', max(bar_len, 0)) * "┫ $b → $d_str"
            end
            push!(lines, line)
        end
    end

    return join(lines, "\n")
end

"""
    plot_ascii(diagram::PersistenceDiagram) -> String

ASCII plot of a persistence diagram. Birth on x-axis, death on y-axis.
"""
function plot_ascii(diagram::PersistenceDiagram)
    isempty(diagram.points) && return "(empty diagram)"

    max_val = 0.0
    for (b, d) in diagram.points
        max_val = max(max_val, isinf(d) ? b + 1.0 : d)
    end

    h = 20
    w = 40
    grid = fill(' ', h, w)

    # Draw diagonal
    for i in 1:min(h, w)
        grid[h - i + 1, i] = '\\'
    end

    # Plot points
    for (b, d) in diagram.points
        if isinf(d)
            # Plot as point on top edge
            bx = clamp(round(Int, b / max_val * (w - 1)) + 1, 1, w)
            grid[1, bx] = '★'
        else
            bx = clamp(round(Int, b / max_val * (w - 1)) + 1, 1, w)
            by = clamp(h - round(Int, d / max_val * (h - 1)), 1, h)
            grid[by, bx] = '●'
        end
    end

    lines = String[]
    push!(lines, "  death ↑")
    for r in 1:h
        row_str = String(grid[r, :])
        push!(lines, "  " * row_str)
    end
    push!(lines, "  " * repeat('─', w))
    push!(lines, "  birth →")

    return join(lines, "\n")
end

import Base: show
show(io::IO, b::PersistenceBarcode) = print(io, "PersistenceBarcode($(num_bars(b)) bars)")
show(io::IO, d::PersistenceDiagram) = print(io, "PersistenceDiagram($(num_points(d)) points)")
