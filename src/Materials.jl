abstract type AbstractMaterial end

struct Material{MPT} <: AbstractMaterial
    ID
    E::MPT
    ν::MPT

    function Material(ID, E::Real, ν::Real)
        material_properties = promote(E, ν)

        return new{eltype(material_properties)}(ID, material_properties...)
    end
end