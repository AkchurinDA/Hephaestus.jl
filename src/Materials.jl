"""
    struct Material

A type representing a material in the finite element model of a structure of interest.
"""
struct Material{T<:Real}
    ID    ::Int
    E     ::T
    ν     ::T
    ρ     ::T
end

Material(ID::Int, E::Real, ν::Real, ρ::Real) = Material(ID, promote(E, ν, ρ)...)