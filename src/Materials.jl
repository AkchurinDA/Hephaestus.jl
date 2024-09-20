"""
    Material

A type that represents a material in the finite element model of a structure.

# Fields
$(FIELDS)
"""
struct Material{MPT<:Real}
    "Unique identifier of the material provided by the user"
    ID  ::Int
    "Young's modulus of the material, ``E``"
    E   ::MPT
    "Poisson's ratio of the material, ``\\nu``"
    ν   ::MPT
    "Density of the material, ``\\rho``"
    ρ   ::MPT
end

Material(ID::Int, E::Real, ν::Real, ρ::Real) = Material(ID, promote(E, ν, ρ)...)