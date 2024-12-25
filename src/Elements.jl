struct Element{NIT<:Real, NJT<:Real, ST<:Real, MT<:Real}
    "Identification tag"
    ID::Int
    "Node (``i``)"
    node_i::Node{NIT}
    "Node (``j``)"
    node_j::Node{NJT}
    "Section"
    section::Section{ST}
    "Material"
    material::Material{MT}
    "Orientation angle, ``\\omega``"
    ω::Real
end