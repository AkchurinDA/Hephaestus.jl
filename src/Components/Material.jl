struct Material{T <: Real}
    "Material ID."
    ID::Integer
    "Elastic modulus."
    E::T
    "Poisson's ratio."
    v::T
    "Shear modulus."
    G::T

    function Material(ID::Integer, E::Real, v::Real)
        # Compute the shear modulus:
        G = E / (2 * (1 + v))

        # Promote the material properties to the same common type:
        MaterialProperties = promote(E, v, G)

        # Construct the material:
        return new{eltype(MaterialProperties)}(ID, MaterialProperties...)
    end
end