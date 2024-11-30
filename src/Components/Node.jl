struct Node{T<:Real} <: AbstractNode
    x::T
    y::T
    z::T

    function Node(x::Real, y::Real, z::Real)
        T = float(promote_type(
            typeof(x), 
            typeof(y), 
            typeof(z)))

        return new{T}(x, y, z)
    end
end