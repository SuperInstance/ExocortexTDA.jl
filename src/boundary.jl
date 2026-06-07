# ============================================================
# boundary.jl — Boundary matrix reduction
# ============================================================

"""
    BoundaryMatrix

Sparse representation of a mod-2 boundary matrix.
Columns correspond to simplices; rows correspond to (dim-1)-simplices.
Entries are 0 or 1 (mod 2).
"""
struct BoundaryMatrix
    # For each column j, the set of row indices with 1s
    columns::Vector{Set{Int}}
    simplex_order::Vector{Simplex}
    index_map::Dict{Simplex, Int}
end

"""
    boundary_matrix(sc::SimplicialComplex)

Build mod-2 boundary matrix from a simplicial complex.
Rows and columns are ordered by (dimension, vertices).
Returns a BoundaryMatrix with:
- columns: sparse column representations
- simplex_order: the ordering of simplices (column order)
- index_map: simplex → column index
"""
function boundary_matrix(sc::SimplicialComplex)
    all_simps = collect(sc.simplices)
    # Sort by (dimension, vertices) for consistent ordering
    sort!(all_simps, by=s -> (dim(s), vertices(s)))

    n = length(all_simps)
    index_map = Dict{Simplex, Int}()
    for (i, s) in enumerate(all_simps)
        index_map[s] = i
    end

    columns = Vector{Set{Int}}(undef, n)
    for j in 1:n
        s = all_simps[j]
        faces = boundary_mod2(s)
        col = Set{Int}()
        for face in faces
            if haskey(index_map, face)
                push!(col, index_map[face])
            end
        end
        columns[j] = col
    end

    return BoundaryMatrix(columns, all_simps, index_map)
end

"""
    boundary_matrix(f::Filtration)

Build boundary matrix from a filtration (uses the filtration ordering).
"""
function boundary_matrix(f::Filtration)
    n = length(f.simplices)
    index_map = Dict{Simplex, Int}()
    for (i, s) in enumerate(f.simplices)
        index_map[s] = i
    end

    columns = Vector{Set{Int}}(undef, n)
    for j in 1:n
        s = f.simplices[j]
        faces = boundary_mod2(s)
        col = Set{Int}()
        for face in faces
            if haskey(index_map, face)
                push!(col, index_map[face])
            end
        end
        columns[j] = col
    end

    return BoundaryMatrix(columns, f.simplices, index_map)
end

"""Lowest nonzero row index in column j. Returns 0 if empty."""
function low(M::BoundaryMatrix, j::Int)
    isempty(M.columns[j]) && return 0
    return maximum(M.columns[j])
end

"""
    reduce_matrix!(M::BoundaryMatrix) -> Vector{Tuple{Int, Int}}

Reduce boundary matrix using iterative standard algorithm (mod 2).
Returns vector of (birth, death) column pairs representing persistence pairs.
Also identifies unpaired columns (essential cycles).

Algorithm (iterative):
```
for j = 1..n:
    while low(j) != 0 and low(j) was seen:
        k = column that has low(j)
        M[:,j] = M[:,j] ⊕ M[:,k]   (mod 2)
    if low(j) != 0:
        mark low(j) as paired with j
```
"""
function reduce_matrix!(M::BoundaryMatrix)
    n = length(M.columns)
    # Track which low value maps to which column
    low_to_col = Dict{Int, Int}()
    pairs = Tuple{Int, Int}[]

    for j in 1:n
        while true
            l = low(M, j)
            l == 0 && break
            if haskey(low_to_col, l)
                # XOR with the column that has the same low
                k = low_to_col[l]
                # Symmetric difference (mod 2 addition)
                M.columns[j] = symmetric_diff(M.columns[j], M.columns[k])
            else
                break
            end
        end
        l = low(M, j)
        if l != 0
            low_to_col[l] = j
            push!(pairs, (l, j))
        end
    end

    return pairs
end

"""Symmetric difference of two sets (mod-2 addition)."""
function symmetric_diff(a::Set{Int}, b::Set{Int})
    result = Set{Int}()
    for x in a
        if x in b
            # Cancel out (mod 2)
        else
            push!(result, x)
        end
    end
    for x in b
        if !(x in a)
            push!(result, x)
        end
    end
    return result
end

"""Get unpaired (essential) column indices."""
function unpaired_columns(M::BoundaryMatrix, pairs::Vector{Tuple{Int, Int}})
    paired = Set{Int}()
    for (b, d) in pairs
        push!(paired, b)
        push!(paired, d)
    end
    unpaired = Int[]
    for j in 1:length(M.columns)
        if !(j in paired)
            push!(unpaired, j)
        end
    end
    return unpaired
end
