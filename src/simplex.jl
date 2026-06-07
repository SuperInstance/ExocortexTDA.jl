# ============================================================
# simplex.jl — Simplicial complex operations
# ============================================================

"""
    Simplex

A simplex represented by a sorted vector of integer vertices.

Multiple dispatch constructors:
- `Simplex(a::Int)` → 0-simplex (vertex)
- `Simplex(a::Int, b::Int)` → 1-simplex (edge)
- `Simplex(vertices::Vector{Int})` → k-simplex
"""
struct Simplex
    verts::Vector{Int}
    function Simplex(verts::Vector{Int})
        v = sort(unique(verts))
        @assert length(v) == length(verts) "Simplex vertices must be unique"
        new(v)
    end
end

Simplex(a::Int) = Simplex([a])
Simplex(a::Int, b::Int) = Simplex(a == b ? [a] : [a, b])

"""Dimension of a simplex."""
dim(s::Simplex) = length(s.verts) - 1

"""Sorted vertex list of a simplex."""
vertices(s::Simplex) = s.verts

import Base: ==
==(a::Simplex, b::Simplex) = a.verts == b.verts
import Base: hash
hash(s::Simplex, h::UInt) = hash(s.verts, h)
import Base: isless
isless(a::Simplex, b::Simplex) = (dim(a), a.verts) < (dim(b), b.verts)
import Base: show
show(io::IO, s::Simplex) = print(io, "Simplex(", join(s.verts, ","), ") [dim=$(dim(s))]")

"""
    boundary(s::Simplex) -> Vector{Tuple{Simplex, Int}}

Compute the boundary of a simplex as a list of (face, sign) pairs.
The sign alternates by the position of the omitted vertex (mod-2 ignored here; sign kept for reference).
"""
function boundary(s::Simplex)
    n = length(s.verts)
    n == 1 && return Tuple{Simplex, Int}[]
    result = Tuple{Simplex, Int}[]
    for i in 1:n
        face_verts = vcat(s.verts[1:i-1], s.verts[i+1:end])
        face = Simplex(face_verts)
        push!(result, (face, (-1)^(i - 1)))
    end
    return result
end

"""
    boundary_mod2(s::Simplex) -> Vector{Simplex}

Compute mod-2 boundary (no signs).
"""
function boundary_mod2(s::Simplex)
    n = length(s.verts)
    n == 1 && return Simplex[]
    result = Simplex[]
    for i in 1:n
        face_verts = vcat(s.verts[1:i-1], s.verts[i+1:end])
        push!(result, Simplex(face_verts))
    end
    return result
end

"""
    SimplicialComplex

A simplicial complex stored as a Set of Simplex objects, plus a sorted index by dimension.
"""
mutable struct SimplicialComplex
    simplices::Set{Simplex}
    by_dim::Dict{Int, Set{Simplex}}

    SimplicialComplex() = new(Set{Simplex}(), Dict{Int, Set{Simplex}}())
end

function SimplicialComplex(simplexes::AbstractVector{Simplex})
    sc = SimplicialComplex()
    for s in simplexes
        add_simplex!(sc, s)
    end
    return sc
end

"""Add a simplex and all its faces to the complex."""
function add_simplex!(sc::SimplicialComplex, s::Simplex)
    push!(sc.simplices, s)
    d = dim(s)
    if !haskey(sc.by_dim, d)
        sc.by_dim[d] = Set{Simplex}()
    end
    push!(sc.by_dim[d], s)
    # Add all faces to ensure closure
    if d > 0
        for face in boundary_mod2(s)
            if !(face in sc.simplices)
                add_simplex!(sc, face)
            end
        end
    end
    return sc
end

"""All simplices in the complex."""
simplices(sc::SimplicialComplex) = collect(sc.simplices)

"""Simplices of a given dimension."""
simplices(sc::SimplicialComplex, d::Int) = haskey(sc.by_dim, d) ? collect(sc.by_dim[d]) : Simplex[]

"""All faces of dimension d (alias for simplices(sc, d))."""
faces(sc::SimplicialComplex, d::Int) = simplices(sc, d)

"""Check if simplex is in the complex."""
has_simplex(sc::SimplicialComplex, s::Simplex) = s in sc.simplices

import Base: length
"""Number of simplices."""
length(sc::SimplicialComplex) = length(sc.simplices)

"""Maximum dimension of any simplex."""
max_dim(sc::SimplicialComplex) = isempty(sc.by_dim) ? -1 : maximum(keys(sc.by_dim))

"""Number of simplices at each dimension (f-vector)."""
f_vector(sc::SimplicialComplex) = [length(simplices(sc, d)) for d in 0:max_dim(sc)]

"""
    euler_characteristic(sc::SimplicialComplex)

Compute the Euler characteristic: Σ (-1)^k * f_k
"""
function euler_characteristic(sc::SimplicialComplex)
    md = max_dim(sc)
    ec = 0
    for d in 0:md
        ec += (-1)^d * length(simplices(sc, d))
    end
    return ec
end
