"""
    struct Node

A type that represents a node in the FE model of a structure.

This type should not be used directly. 
Instead, use the [`add_node!()`](@ref) function to add a node to the model.
To remove a node from the model, use the [`del_node!()`](@ref) function.
"""
mutable struct Node{CT <: Real}
    "Unique identifier of the node provided by the user"
    ID_u            ::Int
    "Unique identifier of the node assigned by the program"
    ID_i            ::Union{Nothing, <:Int}
    "``x``-coordinate"
    x               ::CT
    "``y``-coordinate"
    y               ::CT
    "``z``-coordinate"
    z               ::CT
    "Displacement in the global ``x``-direction"
    u_x             ::Real
    "Displacement in the global ``y``-direction"
    u_y             ::Real
    "Displacement in the global ``z``-direction"
    u_z             ::Real
    "Rotation about the global ``x``-axis"
    θ_x             ::Real
    "Rotation about the global ``y``-axis"
    θ_y             ::Real
    "Rotation about the global ``z``-axis"
    θ_z             ::Real
    "Is the nodal translation fixed in the global ``x``-direction?"
    u_x_supported   ::Bool
    "Is the nodal translation fixed in the global ``y``-direction?"
    u_y_supported   ::Bool
    "Is the nodal translation fixed in the global ``z``-direction?"
    u_z_supported   ::Bool
    "Is the nodal rotation fixed about the global ``x``-axis?"
    θ_x_supported   ::Bool
    "Is the nodal rotation fixed about the global ``y``-axis?"
    θ_y_supported   ::Bool
    "Is the nodal rotation fixed about the global ``z``-axis?"
    θ_z_supported   ::Bool
    "Nodal load in the global ``x``-direction"
    F_x             ::Real
    "Nodal load in the global ``y``-direction"
    F_y             ::Real
    "Nodal load in the global ``z``-direction"
    F_z             ::Real
    "Nodal moment about the global ``x``-axis"
    M_x             ::Real
    "Nodal moment about the global ``y``-axis"
    M_y             ::Real
    "Nodal moment about the global ``z``-axis"
    M_z             ::Real
    u_x_enforced    ::Real
    u_y_enforced    ::Real
    u_z_enforced    ::Real
    θ_x_enforced    ::Real
    θ_y_enforced    ::Real
    θ_z_enforced    ::Real

    function Node(ID::Int, x::CT, y::CT, z::CT) where {CT <: Real}        
        new{CT}(ID, nothing,
            x, y, z, 
            fill(0    , 6)...,
            fill(false, 6)...,
            fill(0    , 6)...,
            fill(0    , 6)...)
    end
end

Node(ID::Int, x::Real, y::Real, z::Real) = Node(ID, promote(x, y, z)...)