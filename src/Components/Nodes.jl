"""
    struct Node

A type representing a node in a finite element model.

$(FIELDS)
"""
struct Node{T<:Real}
    "Identification tag"
    ID::Int
    "``x``-coordinate"
    x::T
    "``y``-coordinate"
    y::T
    "``z``-coordinate"
    z::T
    "Is the node restrained agaist translation in the global ``x``-direction?"
    u_x::Bool
    "Is the node restrained agaist translation in the global ``y``-direction?"
    u_y::Bool
    "Is the node restrained agaist translation in the global ``z``-direction?"
    u_z::Bool
    "Is the node restrained agaist rotation about the global ``x``-axis?"
    θ_x::Bool
    "Is the node restrained agaist rotation about the global ``y``-axis?"
    θ_y::Bool
    "Is the node restrained agaist rotation about the global ``z``-axis?"
    θ_z::Bool

    function Node(ID::Int, 
        x::T1, y::T2, z::T3,
        u_x::Bool, u_y::Bool, u_z::Bool,
        θ_x::Bool, θ_y::Bool, θ_z::Bool) where {
        T1<:Real,
        T2<:Real,
        T3<:Real}
        # Promote the types:
        T = float(promote_type(T1, T2, T3))

        # Construct a new instance:
        return new{T}(ID, x, y, z, u_x, u_y, u_z, θ_x, θ_y, θ_z)
    end
end

