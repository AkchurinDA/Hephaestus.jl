abstract type AbstractNode end

struct Node{CT} <: AbstractNode
    ID
    x::CT
    y::CT
    z::CT

    function Node(ID, x::Real, y::Real, z::Real)
        coordinates = promote(x, y, z)

        return new{eltype(coordinates)}(ID, coordinates...)
    end
end