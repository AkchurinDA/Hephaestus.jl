"""
    getnodaldisplacements(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the displacement vector of a node of interest from the solution cache.
"""
getnodaldisplacements(model::Model, solution::AbstractSolutionCache, ID::Int) = getnodaldisplacements(model, solution.U, ID)
function getnodaldisplacements(model::Model, U::AbstractVector{<:Real}, ID::Int)
    # Find the index of the node in the model:
    index = findfirst(x -> x.ID == ID, model.nodes) 

    # Extract the displacement vector of the node:
    nodedisp = U[(6 * index - 5):(6 * index)]

    # Return the displacement vector:
    return nodedisp
end

"""
    getnodalreactions(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the reaction vector of a node of interest from the solution cache.
"""
getnodalreactions(model::Model, solution::AbstractSolutionCache, ID::Int) = getnodalreactions(model, solution.R, ID)
function getnodalreactions(model::Model, R::AbstractVector{<:Real}, ID::Int)
    # Find the index of the node in the model:
    index = findfirst(x -> x.ID == ID, model.nodes) 

    # Extract the reaction vector of the node:
    nodereactions = R[(6 * index - 5):(6 * index)]

    # Return the reaction vector:
    return nodereactions
end

"""
    getelementdisplacements(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the displacement vector of an element of interest in its local coordinate system from the solution cache.
"""
getelementdisplacements(model::Model, solution::AbstractSolutionCache, ID::Int) = getelementdisplacements(model, solution.U, ID)
function getelementdisplacements(model::Model, U::AbstractVector{<:Real}, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the displacements of the nodes of the element:
    nodeidisplacements = getnodaldisplacements(model, U, element.node_i.ID)
    nodejdisplacements = getnodaldisplacements(model, U, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    elementdisplacements = [nodeidisplacements; nodejdisplacements]

    # Transform the displacement vector of the element to the local coordinate system:
    elementdisplacements = element.Γ * elementdisplacements

    # Return the displacement vector of the element:
    return elementdisplacements
end

"""
    getelementforces(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the internal forces of an element of interest in its local coordinate system from the solution cache.
"""
getelementforces(model::Model, solution::AbstractSolutionCache, ID::Int) = getelementforces(model, solution.U, ID)
function getelementforces(model::Model, U::AbstractVector{<:Real}, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the displacements of the nodes of the element:
    nodeidisp = getnodaldisplacements(model, U, element.node_i.ID)
    nodejdisp = getnodaldisplacements(model, U, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    elementdisplacements = [nodeidisp; nodejdisp]

    # Transform the displacement vector of the element to the local coordinate system:
    elementdisplacements = element.Γ * elementdisplacements

    # Compute the internal forces of the element:
    elementforces = element.k_e_l * elementdisplacements

    # Return the internal forces of the element:
    return elementforces
end

function getelementforces(element::Element, elementstate::ElementState)
    N = elementstate.N
    k_e_l = element.k_e_l
    k_g_l = N * element.k_g_l

    elementdisplacements = [elementstate.node_i_coords; elementstate.node_j_coords]

    elementdisplacements = elementstate.Γ * elementdisplacements
    
    elementforces = (k_e_l + k_g_l) * elementdisplacements

    return elementforces
end