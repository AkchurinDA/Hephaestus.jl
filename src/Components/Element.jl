struct Element{NICT<:Real, NJCT<:Real, MPT<:Real, SPT<:Real} <: AbstractElement
    node_i_ID::Int
    x_i::NICT
    y_i::NICT
    z_i::NICT
    node_j_ID::Int
    x_j::NJCT
    y_j::NJCT
    z_j::NJCT
    material_ID::Int
    E::MPT
    ν::MPT
    ρ::MPT
    section_ID::Int
    A   ::SPT
    I_zz::SPT
    I_yy::SPT
    J   ::SPT

    ω::Real
end