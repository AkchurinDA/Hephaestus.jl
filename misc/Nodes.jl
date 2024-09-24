"""
    struct Node

A type that represents a node in the finite element model of a structure.

This type should never be called directly by the user.

# Fields
$(FIELDS)
"""
struct Node{CT<:Real}
    "Unique identifier of the node provided by the user"
    ID  ::Int
    "``x``-coordinate of the node"
    x   ::CT
    "``y``-coordinate of the node"
    y   ::CT
    "``z``-coordinate of the node"
    z   ::CT
end

Node(ID::Int, x::Real, y::Real, z::Real) = Node(ID, promote(x, y, z)...)