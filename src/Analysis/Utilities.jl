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
    node_i_disp = extract_node_disp(model, U, element.node_i.ID)
    node_j_disp = extract_node_disp(model, U, element.node_j.ID)

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
    node_i_disp = extract_node_disp(model, U, element.node_i.ID)
    node_j_disp = extract_node_disp(model, U, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    element_disp_g = [node_i_disp; node_j_disp]

    # Return the displacement vector of the element:
    return element_disp_g
end