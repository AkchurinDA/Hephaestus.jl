"""
    struct ConcentratedLoad

A type representing a concentrated load in the finite element model.

$(FIELDS)
"""
struct ConcentratedLoad{T<:Real}
    "Identification tag of the node at which the concentrated load is applied"
    ID::Int
    "Concentrated load in the global ``x``-direction"
    F_x::T
    "Concentrated load in the global ``y``-direction"
    F_y::T
    "Concentrated load in the global ``z``-direction"
    F_z::T
    "Moment about the global ``x``-axis"
    M_x::T
    "Moment about the global ``y``-axis"
    M_y::T
    "Moment about the global ``z``-axis"
    M_z::T

    function ConcentratedLoad(ID::Int, 
        F_x::T1, F_y::T2, F_z::T3, 
        M_x::T4, M_y::T5, M_z::T6) where {
        T1<:Real,
        T2<:Real,
        T3<:Real,
        T4<:Real,
        T5<:Real,
        T6<:Real}
        # Promote the types:
        T = float(promote_type(T1, T2, T3, T4, T5, T6))

        # Construct a new instance:
        return new{T}(ID, F_x, F_y, F_z, M_x, M_y, M_z)
    end
end

get_concload_T(::ConcentratedLoad{T}) where {T} = T

"""
    struct DistributedLoad

A type representing a concentrated load in the finite element model.

$(FIELDS)
"""
struct DistributedLoad{T<:Real}
    "Identification tag of the element to which the distributed load is applied"
    ID::Int
    "Distributed load in the local/global ``x``-direction"
    w_x::T
    "Distributed load in the local/global ``y``-direction"
    w_y::T
    "Distributed load in the local/global ``z``-direction"
    w_z::T
    "Fixed-end force vector in the global coordinate system"
    p::Vector{T}

    function DistributedLoad(ID::Int, 
        w_x::T1, w_y::T2, w_z::T3,
        p::AbstractVector{T4}) where {
        T1<:Real,
        T2<:Real,
        T3<:Real,
        T4<:Real}
        # Promote the types:
        T = float(promote_type(T1, T2, T3, T4))
        
        # Construct a new instance:
        return new{T}(ID, w_x, w_y, w_z, p)
    end
end

get_distload_T(::DistributedLoad{T}) where {T} = T