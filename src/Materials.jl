struct Material{MPT<:Real}
    tag ::Int
    E   ::MPT
    ν   ::MPT
    ρ   ::MPT
end

Material(tag::Int, E::Real, ν::Real, ρ::Real) = Material(tag, promote(E, ν, ρ)...)