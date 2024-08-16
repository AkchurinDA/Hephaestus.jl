abstract type AbstractElement end

struct Element{NICT,NJCT,MPT,SPT,EOT,EPT} <: AbstractElement
    ID
    "1st node of the element"
    node_i::Node{NICT}
    "2nd node of the element"
    node_j::Node{NJCT}
    "Material properties of the element"
    material::Material{MPT}
    "Section properties of the element"
    section::Section{SPT}
    "Orientation of the element"
    ω::EOT
    "Length of the element"
    L::EPT

    function Element(ID, node_i::Node{NICT}, node_j::Node{NJCT}, material::Material{MPT}, section::Section{SPT}, ω::EOT) where {NICT,NJCT,MPT,SPT,EOT}
        # Compute the length of the element:
        L = compute_L(
            node_i.x, node_i.y, node_i.z,
            node_j.x, node_j.y, node_j.z)

        # Return the element:
        return new{NICT, NJCT, MPT, SPT, EOT, typeof(L)}(ID, node_i, node_j, material, section, L)
    end
end

function compute_L(
    x_i::NICT, y_i::NICT, z_i::NICT,
    x_j::NJCT, y_j::NJCT, z_j::NJCT) where {NICT,NJCT}
    # Compute the length of the element:
    L = sqrt((x_j - x_i) ^ 2 + (y_j - y_i) ^ 2 + (z_j - z_i) ^ 2)

    # Return the result:
    return L
end

