struct Material{T <: Real}
    "Material ID."
    ID::Integer
    "Elastic modulus."
    E::T

    function Material(ID::Integer, E::Real)
        # Promote the material properties to the same common type:
        material_properties = promote(E)

        # Construct the material:
        return new{eltype(material_properties)}(ID, material_properties...)
    end
end