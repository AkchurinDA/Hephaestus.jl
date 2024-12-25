function solve end

abstract type AbstractAnalysisType end
include("AnalysisTypes/LinearElasticAnalysis.jl")
include("AnalysisTypes/NonlinearElasticAnalysis.jl")

function assemble_K_e(model::Model)
    # Initialize the global elastic stiffness matrix:
    K_e = zeros(6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global elastic stiffness matrix:
    for element in model.elements
        # Extact the transformation matrix of the element:
        Γ = element.Γ

        # Compute the element stiffness matrix in the global coordinate system:
        k_e_g = Γ' * element.k_e_l * Γ

        # Condense the element stiffness matrix if end releases are present:

        # Assemble the element stiffness matrix into the global stiffness matrix:
        idx_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        idx_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds K_e[(6 * idx_i - 5):(6 * idx_i), (6 * idx_i - 5):(6 * idx_i)] += k_e_g[ 1:6,  1:6]
        @inbounds K_e[(6 * idx_i - 5):(6 * idx_i), (6 * idx_j - 5):(6 * idx_j)] += k_e_g[ 1:6, 7:12]
        @inbounds K_e[(6 * idx_j - 5):(6 * idx_j), (6 * idx_i - 5):(6 * idx_i)] += k_e_g[7:12,  1:6]
        @inbounds K_e[(6 * idx_j - 5):(6 * idx_j), (6 * idx_j - 5):(6 * idx_j)] += k_e_g[7:12, 7:12]
    end

    # Return the global elastic stiffness matrix:
    return K_e
end

function assemble_K_g(model::Model, P::Vector{<:Real})
    # Make sure the length of the vector of axial loads is equal to the number of elements:
    @assert length(P) == length(model.elements) "The length of the vector of axial loads is not equal to the number of elements. Something is wrong."

    # Initialize the global geometric stiffness matrix:
    K_g = zeros(6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global geometric stiffness matrix:
    for (element, p) in zip(model.elements, P)
        # Extact the transformation matrix of the element:
        Γ = element.Γ

        # Compute the element stiffness matrix in the global coordinate system:
        k_e_g = p * Γ' * element.k_e_l * Γ

        # Condense the element stiffness matrix if end releases are present:

        # Assemble the element stiffness matrix into the global stiffness matrix:
        idx_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        idx_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds K_g[(6 * idx_i - 5):(6 * idx_i), (6 * idx_i - 5):(6 * idx_i)] += k_e_g[ 1:6,  1:6]
        @inbounds K_g[(6 * idx_i - 5):(6 * idx_i), (6 * idx_j - 5):(6 * idx_j)] += k_e_g[ 1:6, 7:12]
        @inbounds K_g[(6 * idx_j - 5):(6 * idx_j), (6 * idx_i - 5):(6 * idx_i)] += k_e_g[7:12,  1:6]
        @inbounds K_g[(6 * idx_j - 5):(6 * idx_j), (6 * idx_j - 5):(6 * idx_j)] += k_e_g[7:12, 7:12]
    end

    # Return the global geometric stiffness matrix:
    return K_g
end

function assemble_M(model::Model)
    # Initialize the global mass matrix:
    M = zeros(6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global mass matrix:
    for element in model.elements
        # Extact the transformation matrix of the element:
        Γ = element.Γ

        # Compute the element mass matrix in the global coordinate system:
        m_g = Γ' * element.m_l * Γ

        # Condense the element mass matrix if end releases are present:

        # Assemble the element mass matrix into the global mass matrix:
        idx_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        idx_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds M[(6 * idx_i - 5):(6 * idx_i), (6 * idx_i - 5):(6 * idx_i)] += m_g[ 1:6,  1:6]
        @inbounds M[(6 * idx_i - 5):(6 * idx_i), (6 * idx_j - 5):(6 * idx_j)] += m_g[ 1:6, 7:12]
        @inbounds M[(6 * idx_j - 5):(6 * idx_j), (6 * idx_i - 5):(6 * idx_i)] += m_g[7:12,  1:6]
        @inbounds M[(6 * idx_j - 5):(6 * idx_j), (6 * idx_j - 5):(6 * idx_j)] += m_g[7:12, 7:12]
    end
end

function assemble_F(model::Model)
    # Initialize the global load vector:
    F = zeros(6 * length(model.nodes))

    for concload in model.concloads
        # Extract the concentrated load vector in the global coordinate system:
        F_x, F_y, F_z = concload.F_x, concload.F_y, concload.F_z
        M_x, M_y, M_z = concload.M_x, concload.M_y, concload.M_z

        # Assemble the concentrated load vector into the global load vector:
        idx = findfirst(x -> x.ID == concload.ID, model.nodes)
        @inbounds F[6 * idx - 5] += F_x
        @inbounds F[6 * idx - 4] += F_y
        @inbounds F[6 * idx - 3] += F_z
        @inbounds F[6 * idx - 2] += M_x
        @inbounds F[6 * idx - 1] += M_y
        @inbounds F[6 * idx    ] += M_z
    end

    for distload in model.distloads
        # Compute the distributed load vector in the global coordinate system:

        # Assemble the distributed load vector into the global load vector:
    end

    # Return the global load vector:
    return F
end

function getpartitionindices(model::Model)
    # Initialize the partition indices:
    indices_f = fill(false, 6 * length(model.nodes))
    indices_s = fill(false, 6 * length(model.nodes))

    # Assemble the partition indices:
    for (i, node) in enumerate(model.nodes)
        u_x, u_y, u_z = node.u_x, node.u_y, node.u_z
        θ_x, θ_y, θ_z = node.θ_x, node.θ_y, node.θ_z

        # Assemble the partition indices for the free DOFs:
        indices_f[6 * i - 5] = !u_x
        indices_f[6 * i - 4] = !u_y
        indices_f[6 * i - 3] = !u_z
        indices_f[6 * i - 2] = !θ_x
        indices_f[6 * i - 1] = !θ_y
        indices_f[6 * i    ] = !θ_z

        # Assemble the partition indices for the supported DOFs:
        indices_s[6 * i - 5] = u_x
        indices_s[6 * i - 4] = u_y
        indices_s[6 * i - 3] = u_z
        indices_s[6 * i - 2] = θ_x
        indices_s[6 * i - 1] = θ_y
        indices_s[6 * i    ] = θ_z
    end

    # Return the partition indices:
    return indices_f, indices_s
end

function extract_node_disp(model::Model, U::Vector{<:Real}, ID::Int)
    # Find the index of the node in the model:
    index = findfirst(x -> x.ID == ID, model.nodes) 

    # Extract the displacement vector of the node:
    node_disp = U[(6 * index - 5):(6 * index)]

    # Return the displacement vector:
    return node_disp
end

function extract_element_disp_l(model::Model, U::Vector{<:Real}, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the displacements of the nodes of the element:
    node_i_disp = extract_node_disp(U, element.node_i.ID)
    node_j_disp = extract_node_disp(U, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    element_disp_g = [node_i_disp; node_j_disp]

    # Transform the displacement vector of the element to the local coordinate system:
    element_disp_l = element.Γ * element_disp_g

    # Return the displacement vector of the element:
    return element_disp_l
end

function extract_element_disp_g(model::Model, U::Vector{<:Real}, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the displacements of the nodes of the element:
    node_i_disp = extract_node_disp(U, element.node_i.ID)
    node_j_disp = extract_node_disp(U, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    element_disp_g = [node_i_disp; node_j_disp]

    # Return the displacement vector of the element:
    return element_disp_g
end