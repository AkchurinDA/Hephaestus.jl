abstract type AbstractNode end

"""
    struct Node

A type representing a node in a finite element model.

$(FIELDS)
"""
struct Node{NCT<:Real} <: AbstractNode
    "``x``-coordinate"
    x::NCT
    "``y``-coordinate"
    y::NCT
    "``z``-coordinate"
    z::NCT
end

Node(x::Real, y::Real, z::Real) = Node(promote(x, y, z)...)