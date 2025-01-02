"""
    getnodaldisplacements(model::Model, ID::Int)

Extracts the displacement vector of a node of interest.
"""
function getnodaldisplacements(model::Model, ID::Int)
    # Find the node in the model:
    node = model.nodes[findfirst(x -> x.ID == ID, model.nodes)]

    # Extract the displacement vector of the node:
    δu = [node.state.u_x, node.state.u_y, node.state.u_z, node.state.θ_x, node.state.θ_y, node.state.θ_z]

    # Return the displacement vector:
    return δu
end

function getnodaldisplacements(model::Model, δU::AbstractVector{<:Real}, ID::Int)
    # Find the index of the node in the model:
    index = findfirst(x -> x.ID == ID, model.nodes)

    # Extract the displacement vector of the node:
    δu = δU[(6 * index - 5):(6 * index)]

    # Return the displacement vector:
    return δu
end

"""
    getnodalreactions(model::Model, ID::Int)

Extracts the reaction vector of a node of interest from the solution cache.
"""
function getnodalreactions(model::Model, ID::Int)
    # Find the node in the model:
    node = model.nodes[findfirst(x -> x.ID == ID, model.nodes)]

    # Extract the displacement vector of the node:
    δr = [node.state.F_r_x, node.state.F_r_y, node.state.F_r_z, node.state.M_r_x, node.state.M_r_y, node.state.M_r_z]

    # Return the displacement vector:
    return δr
end

function getnodalreactions(model::Model, δR::AbstractVector{<:Real}, ID::Int)
    # Find the index of the node in the model:
    index = findfirst(x -> x.ID == ID, model.nodes)

    # Extract the reaction vector of the node:
    δr = δR[(6 * index - 5):(6 * index)]

    # Return the reaction vector:
    return δr
end

"""
    getelementdisplacements(model::Model, ID::Int)

Extracts the displacement vector of an element of interest in its local coordinate system.
"""
function getelementdisplacements(model::Model, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the displacements of the nodes of the element:
    δu_i_g = getnodaldisplacements(model, element.node_i.ID)
    δu_j_g = getnodaldisplacements(model, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    δu_g = [δu_i_g; δu_j_g]

    # Transform the displacement vector of the element to the local coordinate system:
    δu_l = element.state.Γ * δu_g

    # Return the displacement vector of the element:
    return δu_l
end

"""
    getelementforces(model::Model, ID::Int)

Extracts the internal force vector of an element of interest in its local coordinate system.
"""
function getelementforces(model::Model, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the internal forces of the element:
    q = element.state.q

    # Return the internal forces of the element:
    return q
end
