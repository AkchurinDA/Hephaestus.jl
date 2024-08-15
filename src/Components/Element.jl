struct Element
    ID::Integer
    node_i::Node
    node_j::Node
    material::Material
    section::Section
    L
    k_e_l::SparseArrays.SparseMatrixCSC
    # k_g_l::SparseArrays.SparseMatrixCSC
    # K_e_l::SparseArrays.SparseMatrixCSC
    # K_g_l::SparseArrays.SparseMatrixCSC

    function Element(ID::Integer, node_i::Node, node_j::Node, material::Material, section::Section)
        # Compute the length of the element:
        L = _compute_L(node_i, node_j)

        # Compute the element elastic stiffness matrix in the local coordinate system:
        k_e_l = _compute_k_e_l(material, section, L)

        # Compute the element geometric stiffness matrix in the local coordinate system:

        # Compute the transformation matrix:
        
        # Compute the element elastic stiffness matrix in the global coordinate system:

        # Compute the element geometric stiffness matrix in the global coordinate system:

        # Compute the element stiffness matrix in the local coordinate system:
        return new(ID, node_i, node_j, material, section, L, k_e_l)
    end
end

function _compute_L(node_i::Node, node_j::Node)
    ΔX = node_j.X - node_i.X
    ΔY = node_j.Y - node_i.Y
    ΔZ = node_j.Z - node_i.Z
    L = sqrt(ΔX^2 + ΔY^2 + ΔZ^2)

    return L
end

function _compute_k_e_l(material, section, L)
    # Extract the material properties:
    E = material.E
    G = material.G

    # Extract the section properties:
    A = section.A
    I_zz = section.I_zz
    I_yy = section.I_yy
    J = section.J

    # Compute the element stiffness matrix in the local coordinate system:
    k_e_l = SparseArrays.spzeros(12, 12)
    @inline k_e_l[[1, 7], [1, 7]] = (E * A / L) * [+1 -1; -1 +1]
    @inline k_e_l[[2, 6, 8, 12], [2, 6, 8, 12]] = (E * I_zz / L ^ 3) * [+12 +6 * L -12 +6 * L; +6 * L +4 * L ^ 2 -6 * L +2 * L ^ 2; -12 -6 * L +12 -6 * L; +6 * L +2 * L ^ 2 -6 * L +4 * L ^ 2]
    @inline k_e_l[[3, 5, 9, 11], [3, 5, 9, 11]] = (E * I_yy / L ^ 3) * [+12 -6 * L -12 -6 * L; -6 * L +4 * L ^ 2 +6 * L +2 * L ^ 2; -12 +6 * L +12 +6 * L; -6 * L +2 * L ^ 2 +6 * L +4 * L ^ 2]
    @inline k_e_l[[4, 10], [4, 10]] = (G * J / L) * [+1 -1; -1 +1]

    # Return the element stiffness matrix in the local coordinate system:
    return k_e_l
end

function _compute_k_g_l()
    # Compute the element stiffness matrix in the local coordinate system:
    k_g_l = SparseArrays.spzeros(12, 12)

    # Return the element stiffness matrix in the local coordinate system:
    return k_g_l
end

function _compute_T(node_i, node_j, L)
    # Extract the nodal coordinates:
    x_i = node_i.X
    y_i = node_i.Y
    z_i = node_i.Z
    x_j = node_j.X
    y_j = node_j.Y
    z_j = node_j.Z

    # Compute the directional cosines:
    γ_11 = (x_j - x_i) / L
    γ_12 = (y_j - y_i) / L
    γ_13 = (z_j - z_i) / L
    # TODO

    # Return the transformation matrix:
    return T
end