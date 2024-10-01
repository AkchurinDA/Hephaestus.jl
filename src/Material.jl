"""
    struct Material

A type representing a material in the model of a structure of interest.

$(FIELDS)
"""
struct Material{T<:Real}
    "Unique identifier"
    ID    ::Int
    "Young's modulus, ``E``"
    E     ::T
    "Poisson's ratio, ``\nu``"
    ν     ::T
    "Density, ``\\rho``"
    ρ     ::T
end

Material(ID::Int, E::Real, ν::Real, ρ::Real) = Material(ID, promote(E, ν, ρ)...)

