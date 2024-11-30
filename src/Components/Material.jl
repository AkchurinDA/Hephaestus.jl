struct Material{T<:Real} <: AbstractMaterial
    E::T
    ν::T
    ρ::T

    function Material(E::Real, ν::Real, ρ::Real)
        T = float(promote_type(
            typeof(E), 
            typeof(ν), 
            typeof(ρ)))

        return new{T}(E, ν, ρ)
    end
end