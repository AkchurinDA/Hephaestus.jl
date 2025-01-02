mutable struct NodeState
    "Identification tag"
    ID::Int
    "Translation along the global ``x``-axis"
    u_x::Real
    "Translation along the global ``y``-axis"
    u_y::Real
    "Translation along the global ``z``-axis"
    u_z::Real
    "Rotation about the global ``x``-axis"
    θ_x::Real
    "Rotation about the global ``y``-axis"
    θ_y::Real
    "Rotation about the global ``z``-axis"
    θ_z::Real
    "Reaction force along the global ``x``-axis"
    F_r_x::Real
    "Reaction force along the global ``y``-axis"
    F_r_y::Real
    "Reaction force along the global ``z``-axis"
    F_r_z::Real
    "Reaction moment about the global ``x``-axis"
    M_r_x::Real
    "Reaction moment about the global ``y``-axis"
    M_r_y::Real
    "Reaction moment about the global ``z``-axis"
    M_r_z::Real
    "Has the state been modified?"
    modified::Bool

    # Constructor:
    NodeState() = new()
end

"""
    struct Node

A type representing a node in a finite element model.

$(FIELDS)
"""
struct Node{T <: Real}
    "Identification tag"
    ID ::Int
    "``x``-coordinate"
    x  ::T
    "``y``-coordinate"
    y  ::T
    "``z``-coordinate"
    z  ::T
    "Is the node restrained against translation along the global ``x``-axis?"
    u_x::Bool
    "Is the node restrained against translation along the global ``y``-axis?"
    u_y::Bool
    "Is the node restrained against translation along the global ``z``-axis?"
    u_z::Bool
    "Is the node restrained against rotation about the global ``x``-axis?"
    θ_x::Bool
    "Is the node restrained against rotation about the global ``y``-axis?"
    θ_y::Bool
    "Is the node restrained against rotation about the global ``z``-axis?"
    θ_z::Bool
    "Current state"
    state::NodeState

    # Constructor:
    function Node(ID::Int,
        x::T1, y::T2, z::T3,
        u_x::Bool, u_y::Bool, u_z::Bool,
        θ_x::Bool, θ_y::Bool, θ_z::Bool,
        state::NodeState) where {
        T1 <: Real,
        T2 <: Real,
        T3 <: Real}
        # Promote the types:
        T = float(promote_type(T1, T2, T3))

        # Construct a new instance:
        return new{T}(ID, x, y, z, u_x, u_y, u_z, θ_x, θ_y, θ_z, state)
    end
end

function initstate!(node::Node)::Node
    # Initialize the state of the node:
    node.state.ID = node.ID
    node.state.u_x = 0
    node.state.u_y = 0
    node.state.u_z = 0
    node.state.θ_x = 0
    node.state.θ_y = 0
    node.state.θ_z = 0
    node.state.F_r_x = 0
    node.state.F_r_y = 0
    node.state.F_r_z = 0
    node.state.M_r_x = 0
    node.state.M_r_y = 0
    node.state.M_r_z = 0
    node.state.modified = false

    # Return the updated node:
    return node
end

function updatestate!(node::Node, δu::AbstractVector{T}, δr::AbstractVector{T})::Node where {T <: Real}
    # Update the nodal coordinates:
    node.state.u_x     += δu[1]
    node.state.u_y     += δu[2]
    node.state.u_z     += δu[3]
    node.state.θ_x     += δu[4]
    node.state.θ_y     += δu[5]
    node.state.θ_z     += δu[6]
    node.state.F_r_x   += δr[1]
    node.state.F_r_y   += δr[2]
    node.state.F_r_z   += δr[3]
    node.state.M_r_x   += δr[4]
    node.state.M_r_y   += δr[5]
    node.state.M_r_z   += δr[6]
    node.state.modified = true

    # Return the updated node:
    return node
end

gettype(::Node{T}) where {T} = T
