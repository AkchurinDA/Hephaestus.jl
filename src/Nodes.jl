struct Node{CT<:Real}
    tag ::Int
    x   ::CT
    y   ::CT
    z   ::CT
end

Node(tag::Int, x::Real, y::Real, z::Real) = Node(tag, promote(x, y, z)...)