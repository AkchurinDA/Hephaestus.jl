"""
    struct Node

A type representing a node in the finite element model of a structure of interest.
"""
struct Node{T<:Real}
    ID    ::Int
    x     ::T
    y     ::T
    z     ::T
end

Node(ID::Int, x::Real, y::Real, z::Real) = Node(ID, promote(x, y, z)...)