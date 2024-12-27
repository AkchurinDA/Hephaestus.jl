"""
    getnodedisp(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the displacement vector of a node of interest from the solution cache.
"""
getnodedisp(model::Model, solution::AbstractSolutionCache, ID::Int) = getnodedisp(model, solution.U, ID)
function getnodedisp(model::Model, U::AbstractVector{<:Real}, ID::Int)
    # Find the index of the node in the model:
    index = findfirst(x -> x.ID == ID, model.nodes) 

    # Extract the displacement vector of the node:
    nodedisp = U[(6 * index - 5):(6 * index)]

    # Return the displacement vector:
    return nodedisp
end

"""
    getnodereactions(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the reaction vector of a node of interest from the solution cache.
"""
getnodereactions(model::Model, solution::AbstractSolutionCache, ID::Int) = getnodereactions(model, solution.R, ID)
function getnodereactions(model::Model, R::AbstractVector{<:Real}, ID::Int)
    # Find the index of the node in the model:
    index = findfirst(x -> x.ID == ID, model.nodes) 

    # Extract the reaction vector of the node:
    nodereactions = R[(6 * index - 5):(6 * index)]

    # Return the reaction vector:
    return nodereactions
end

"""
    getelementdisp_l(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the displacement vector of an element of interest in the local coordinate system from the solution cache.
"""
getelementdisp_l(model::Model, solution::AbstractSolutionCache, ID::Int) = getelementdisp_l(model, solution.U, ID)
function getelementdisp_l(model::Model, U::AbstractVector{<:Real}, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the displacements of the nodes of the element:
    nodeidisp = getnodedisp(model, U, element.node_i.ID)
    nodejdisp = getnodedisp(model, U, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    elementdisp_g = [nodeidisp; nodejdisp]

    # Transform the displacement vector of the element to the local coordinate system:
    elementdisp_l = element.Γ * elementdisp_g

    # Return the displacement vector of the element:
    return elementdisp_l
end

"""
    getelementdisp_g(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the displacement vector of an element of interest in the global coordinate system from the solution cache.
"""
getelementdisp_g(model::Model, solution::AbstractSolutionCache, ID::Int) = getelementdisp_g(model, solution.U, ID)
function getelementdisp_g(model::Model, U::AbstractVector{<:Real}, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the displacements of the nodes of the element:
    nodeidisp = getnodedisp(model, U, element.node_i.ID)
    nodejdisp = getnodedisp(model, U, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    elementdisp_g = [nodeidisp; nodejdisp]

    # Return the displacement vector of the element:
    return elementdisp_g
end

"""
    getelementforces_l(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the internal forces of an element of interest in the local coordinate system from the solution cache.
"""
getelementforces_l(model::Model, solution::AbstractSolutionCache, ID::Int) = getelementforces_l(model, solution.U, ID)
function getelementforces_l(model::Model, U::AbstractVector{<:Real}, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the displacements of the nodes of the element:
    nodeidisp = getnodedisp(model, U, element.node_i.ID)
    nodejdisp = getnodedisp(model, U, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    elementdisp_g = [nodeidisp; nodejdisp]

    # Transform the displacement vector of the element to the local coordinate system:
    elementdisp_l = element.Γ * elementdisp_g

    # Compute the internal forces of the element:
    elementforces_l = element.k_e_l * elementdisp_l

    # Return the internal forces of the element:
    return elementforces_l
end

"""
    getelementforces_g(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the internal forces of an element of interest in the global coordinate system from the solution cache.
"""
getelementforces_g(model::Model, solution::AbstractSolutionCache, ID::Int) = getelementforces_g(model, solution.U, ID)
function getelementforces_g(model::Model, U::AbstractVector{<:Real}, ID::Int)
    # Find the element in the model:
    element = model.elements[findfirst(x -> x.ID == ID, model.elements)]

    # Extract the displacements of the nodes of the element:
    nodeidisp = getnodedisp(model, U, element.node_i.ID)
    nodejdisp = getnodedisp(model, U, element.node_j.ID)

    # Combine the displacements of the nodes into the displacement vector of the element:
    elementdisp_g = [nodeidisp; nodejdisp]

    # Combine the displacements of the nodes into the displacement vector of the element:
    elementdisp_g = [nodeidisp; nodejdisp]

    # Transform the displacement vector of the element to the local coordinate system:
    elementdisp_l = element.Γ * elementdisp_g

    # Compute the internal forces of the element:
    elementforces_l = element.k_e_l * elementdisp_l
    
    # Transform the internal forces of the element to the global coordinate system:
    elementforces_g = element.Γ * elementforces_l

    # Return the internal forces of the element:
    return elementforces_g
end

"""
    getelementaxialload(model::Model, solution::AbstractSolutionCache, ID::Int)

Extracts the interal axial load of an element of interest from the solution cache.
"""
function getelementaxialload(model::Model, solution::AbstractSolutionCache, ID::Int)
    # Extract the internal forces of the element:
    elementforces_l = getelementforces_l(model, solution, ID)

    # Extract the axial load of the element:
    elementaxialload = elementforces_l[7]

    # Return the axial load of the element:
    return elementaxialload
end