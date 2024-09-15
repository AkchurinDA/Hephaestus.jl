"""
    struct Material

A type that represents a material in the FE model of a structure.

This type should not be used directly. Instead, use the [`add_material!`](@ref) function to add a material to the model.
To remove a material from the model, use the [`del_material!`](@ref) function.
"""
struct Material{MPT <: Real}
    ID  ::Int
    E   ::MPT
    ν   ::MPT
    ρ   ::MPT
end

Material(ID::Int, E::Real, ν::Real, ρ::Real) = Material(ID, promote(E, ν, ρ)...)