struct Node{T <: Real}
    ID::Integer
    X::T
    Y::T
    Z::T

    function Node(ID::Integer, X::Real, Y::Real, Z::Real)
        coordinates = promote(X, Y, Z)

        return new{eltype(coordinates)}(ID, coordinates...)
    end
end