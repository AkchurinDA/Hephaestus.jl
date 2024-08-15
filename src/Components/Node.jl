struct Node{T <: Real}
    ID::Integer
    x::T
    y::T

    function Node(ID::Integer, x::Real, y::Real)
        coordinates = promote(x, y)

        return new{eltype(coordinates)}(ID, coordinates...)
    end
end