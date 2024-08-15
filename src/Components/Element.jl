struct Element
    ID::Integer
    node_i::Node
    node_j::Node
    material::Material
    section::Section
    L::Real
    θ::Real
    k_e_l::SparseArrays.SparseMatrixCSC
    k_g_l::SparseArrays.SparseMatrixCSC
    T::SparseArrays.SparseMatrixCSC
    # K_e_l::SparseArrays.SparseMatrixCSC
    # K_g_l::SparseArrays.SparseMatrixCSC

    function Element(ID::Integer, node_i::Node, node_j::Node, material::Material, section::Section)
        # Compute the length of the element:
        L = _compute_L(node_i, node_j)

        # Compute the orientation of the element:
        θ = _compute_θ(node_i, node_j)

        # Compute the element elastic stiffness matrix in the local coordinate system:
        k_e_l = _compute_k_e_l(material, section, L)

        # Compute the element geometric stiffness matrix in the local coordinate system:
        k_g_l = _compute_k_g_l()

        # Compute the transformation matrix:
        T = _compute_T(θ)
        
        # Compute the element elastic stiffness matrix in the global coordinate system:

        # Compute the element geometric stiffness matrix in the global coordinate system:

        # Compute the element stiffness matrix in the local coordinate system:
        return new(ID, node_i, node_j, material, section, L, θ, k_e_l, k_g_l, T)
    end
end

function _compute_L(node_i::Node, node_j::Node)::Real
    x_i = node_i.x
    y_i = node_i.y
    x_j = node_j.x
    y_j = node_j.y

    Δx = x_j - x_i
    Δy = y_j - y_i

    L = sqrt(Δx ^ 2 + Δy ^ 2)

    return L
end

function _compute_θ(node_i::Node, node_j::Node)::Real
    x_i = node_i.x
    y_i = node_i.y
    x_j = node_j.x
    y_j = node_j.y

    Δx = x_j - x_i
    Δy = y_j - y_i

    θ = atan(Δy, Δx)

    return θ
end

function _compute_k_e_l(material, section, L)::SparseArrays.SparseMatrixCSC
    # Extract the material properties:
    E = material.E

    # Extract the section properties:
    A = section.A
    I = section.I

    # Compute the element stiffness matrix in the local coordinate system:
    k_e_l = SparseArrays.spzeros(6, 6)
    @inline k_e_l[[1, 4], [1, 4]] = (E * A / L) * [+1 -1; -1 +1]
    @inline k_e_l[[2, 3, 5, 6], [2, 3, 5, 6]] = (E * I / L ^ 3) * [+12 +6 * L -12 +6 * L; +6 * L +4 * L ^ 2 -6 * L +2 * L ^ 2; -12 -6 * L +12 -6 * L; +6 * L +2 * L ^ 2 -6 * L +4 * L ^ 2]

    # Return the element stiffness matrix in the local coordinate system:
    return k_e_l
end

function _compute_k_g_l()::SparseArrays.SparseMatrixCSC
    # Compute the element stiffness matrix in the local coordinate system:
    k_g_l = SparseArrays.spzeros(6, 6)

    # Return the element stiffness matrix in the local coordinate system:
    return k_g_l
end

function _compute_T(θ)::SparseArrays.SparseMatrixCSC
    # Compute the transformation matrix:
    T = SparseArrays.spzeros(6, 6)

    # Return the transformation matrix:
    return T
end