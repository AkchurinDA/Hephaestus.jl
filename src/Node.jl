"""
    struct Node

A type representing a node in the model of a structure of interest.

$(FIELDS)
"""
struct Node{T<:Real}
    "Unique identifier"
    ID    ::Int
    "``x``-coordinate"
    x     ::T
    "``y``-coordinate"
    y     ::T
    "``z``-coordinate"
    z     ::T
end

Node(ID::Int, x::Real, y::Real, z::Real) = Node(ID, promote(x, y, z)...)