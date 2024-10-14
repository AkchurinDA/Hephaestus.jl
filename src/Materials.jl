abstract type AbstractMaterial end

"""
    struct Material

A type representing a material in a finite element model.

$(FIELDS)
"""
struct Material{MPT<:Real} <: AbstractMaterial
    "Young's modulus"
    E::MPT
    "Poisson's ratio"
    ν::MPT
    "Density"
    ρ::MPT
end

Material(E::Real, ν::Real, ρ::Real) = Material(promote(E, ν, ρ)...)