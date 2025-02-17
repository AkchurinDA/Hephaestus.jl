"""
    struct Material

A type representing a material in a finite element model.

$(FIELDS)
"""
struct Material{T <: Real}
    "Identification tag"
    ID::Int
    "Young's modulus, ``E``"
    E ::T
    "Poisson's ratio, ``\\nu``"
    ν ::T
    "Density, ``\\rho``"
    ρ ::T

    function Material(ID::Int, 
        E::T1, ν::T2, ρ::T3) where {
        T1 <: Real, 
        T2 <: Real, 
        T3 <: Real}
        # Promote the types:
        T = float(promote_type(T1, T2, T3))

        # Construct a new instance:
        return new{T}(ID, E, ν, ρ)
    end
end

gettype(::Material{T}) where {T} = T