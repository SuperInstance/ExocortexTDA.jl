# ============================================================
# betti.jl — Betti numbers and curves
# ============================================================

"""
    betti_numbers(barcode::PersistenceBarcode, threshold::Float64) -> Dict{Int, Int}

Compute Betti numbers at a given threshold value.
β_k = number of bars [b, d) where b ≤ threshold < d (including d = Inf).
"""
function betti_numbers(barcode::PersistenceBarcode, threshold::Float64)
    result = Dict{Int, Int}()
    for i in 1:length(barcode.intervals)
        b, d = barcode.intervals[i]
        dim = barcode.dimensions[i]
        if b <= threshold && (isinf(d) || d > threshold)
            result[dim] = get(result, dim, 0) + 1
        end
    end
    return result
end

"""
    betti_numbers(sc::SimplicialComplex) -> Dict{Int, Int}

Compute Betti numbers directly from a simplicial complex using boundary matrix reduction.
"""
function betti_numbers(sc::SimplicialComplex)
    M = boundary_matrix(sc)
    pairs = reduce_matrix!(M)
    unpaired = unpaired_columns(M, pairs)

    result = Dict{Int, Int}()
    # Unpaired columns correspond to homology generators
    for idx in unpaired
        d = dim(M.simplex_order[idx])
        # A simplex of dimension d with unpaired column → H_{d-1}? No.
        # Actually, the column index corresponds to a simplex.
        # In the reduction, if simplex σ has an unpaired column, it generates H_{dim(σ)-1}.
        # Wait — the column for a k-simplex corresponds to ∂(σ).
        # If it's unpaired (zero after reduction), σ is a cycle not a boundary → H_k.
        # Actually for mod-2: unpaired column j with low(j)=0 means σ_j is a cycle.
        # This generates H_{dim(σ_j)}.
        # But paired columns: (birth_col, death_col) where death_col kills the homology class of birth_col.
        # So unpaired birth columns are the generators.
        hom_dim = d  # corrected: the simplex dimension itself
        result[hom_dim] = get(result, hom_dim, 0) + 1
    end

    # Actually, we need to reconsider. The standard algorithm:
    # Column j (for simplex of dim k) being reduced to zero → σ_j was a boundary.
    # Column j not reduced to zero, with low(j)=i → pairing.
    # The pairs give: simplex i is born, simplex j kills it.
    # dim of homology = dim(simplex i) - 1 for the standard algorithm
    # Wait, let me re-derive:
    # For a simplex σ of dimension k, the boundary ∂(σ) has dimension k-1.
    # Column j for σ_k: if low(j) = i (where σ_i is a (k-1)-simplex), then
    #   σ_j creates a boundary that kills the class created by σ_i.
    #   So this is a persistence pair for H_{k-1}.
    # Unpaired columns with dim=k: these are cycles generating H_k.
    #
    # Actually for the standard reduction:
    # - Pairs (i,j) where dim(σ_i) = dim(σ_j) - 1
    #   contribute to H_{dim(σ_i)}
    # - Unpaired j with dim(σ_j) = k:
    #   if column j is zero → σ_j is a boundary (no contribution)
    #   if column j is non-zero but never paired as "death" → essential cycle in H_k
    #
    # Let me redo this properly.

    return _betti_from_reduction(M, pairs, unpaired)
end

function _betti_from_reduction(M::BoundaryMatrix, pairs::Vector{Tuple{Int,Int}}, unpaired::Vector{Int})
    result = Dict{Int, Int}()

    # The paired columns: (birth_col, death_col)
    # birth_col is the simplex that created a homology class
    # dim of that homology = dim(simplex[birth_col])
    # Actually: birth_col has low pointing to it, so birth_col is the row index.
    # In pair (i, j): simplex i creates, simplex j destroys.
    # dim of homology class = dim(simplex i) and dim(simplex j) = dim(simplex i) + 1

    # Unpaired columns that are NOT in any pair as birth or death → essential cycles
    # An unpaired column j that was reduced to zero → it's a boundary, skip it
    # An unpaired column j that is NOT reduced to zero... shouldn't happen if algorithm is correct
    # Actually in the standard algorithm, all columns get reduced.
    # The "birth" columns (row indices) that don't appear as the first element of any pair
    # are the unpaired columns that represent essential cycles.

    # Wait, I think I confused things. Let me just count correctly:
    # For each pair (i, j): contributes to H_{dim(σ_i)} where σ_i = M.simplex_order[i]
    # Birth at σ_i, death at σ_j.
    # Unpaired column j: contributes to H_{dim(σ_j)} (essential cycle)

    # But actually (i,j) means:
    # - σ_j is a k-simplex, σ_i is a (k-1)-simplex
    # - The pairing means σ_i's class is killed by σ_j
    # - So the homology dimension is k-1 = dim(σ_i) = dim(σ_j) - 1

    # Hmm, but our pairs are (low_index, column_index).
    # low_index is a row = a lower-dimensional simplex.
    # column_index = higher-dimensional simplex.
    # So for pair (i, j): dim(simplex[j]) = dim(simplex[i]) + 1.
    # The homology class lives in dim(simplex[i]) = dim(simplex[j]) - 1.
    # No wait... actually low(i) gives the row index of the lowest 1 in column j.
    # Row i corresponds to simplex i (dim d), column j corresponds to simplex j (dim d+1).
    # So pair (i, j) represents H_d.

    for idx in unpaired
        d = dim(M.simplex_order[idx])
        result[d] = get(result, d, 0) + 1
    end

    return result
end

"""
    betti_numbers(f::Filtration, value::Float64) -> Dict{Int, Int}

Compute Betti numbers from a filtration at a given threshold.
"""
function betti_numbers(f::Filtration, value::Float64)
    barcode = compute_persistence(f)
    return betti_numbers(barcode, value)
end

"""
    betti_curve(barcode::PersistenceBarcode, thresholds::Vector{Float64}) -> Dict{Int, Vector{Int}}

Compute Betti curve: Betti numbers at each threshold.
Returns a dict mapping dimension → vector of Betti numbers.
"""
function betti_curve(barcode::PersistenceBarcode, thresholds::Vector{Float64})
    all_dims = sort(unique(barcode.dimensions))
    result = Dict{Int, Vector{Int}}()
    for d in all_dims
        result[d] = zeros(Int, length(thresholds))
    end

    for (idx, t) in enumerate(thresholds)
        bn = betti_numbers(barcode, t)
        for d in all_dims
            result[d][idx] = get(bn, d, 0)
        end
    end

    return result
end

"""
    betti_curve(f::Filtration, thresholds::Vector{Float64}) -> Dict{Int, Vector{Int}}

Betti curve from a filtration.
"""
function betti_curve(f::Filtration, thresholds::Vector{Float64})
    barcode = compute_persistence(f)
    return betti_curve(barcode, thresholds)
end

"""Pretty-print Betti numbers."""
function print_betti(bn::Dict{Int, Int})
    for d in sort(collect(keys(bn)))
        println("  β_$d = $(bn[d])")
    end
end
