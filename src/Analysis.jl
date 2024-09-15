function solve(model::Model, type::Symbol)
    # Reset nodal displacements and rotations:
    for node in model.nodes
        node.u_x = 0
        node.u_y = 0
        node.u_z = 0
        node.θ_x = 0
        node.θ_y = 0
        node.θ_z = 0
    end

    internal_IDs_nodes, internal_IDs_elements = _assign_internal_ID(model)
end

function _assign_internal_ID(model::Model)
    # Assign internal IDs to all nodes:
    # User ID => Internal ID
    internal_IDs_nodes = Dict{Int, Int}()
    for (i, node) in enumerate(model.nodes)
        internal_IDs_nodes[node.ID] = i
    end

    # Assign internal IDs to all elements:
    # User ID => Internal ID
    internal_IDs_elements = Dict{Int, Int}()
    for (i, element) in enumerate(model.elements)
        internal_IDs_elements[element.ID] = i
    end

    # Return the internal IDs:
    return internal_IDs_nodes, internal_IDs_elements
end
