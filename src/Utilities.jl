function Base.show(io::IO, model::Model)
    if isempty(model.nodes) && isempty(model.materials) && isempty(model.sections) && isempty(model.elements) && isempty(model.supports) && isempty(model.conc_loads) && isempty(model.dist_loads)
        println(io, styled"{yellow:Empty model}")
    else
        println(io, styled"{yellow:Model with}")
        isempty(model.nodes    ) || println(io, styled"{emphasis:$(length(model.nodes    )) \t Nodes    }")
        isempty(model.materials) || println(io, styled"{emphasis:$(length(model.materials)) \t Materials}")
        isempty(model.sections ) || println(io, styled"{emphasis:$(length(model.sections )) \t Sections }")
        isempty(model.elements ) || println(io, styled"{emphasis:$(length(model.elements )) \t Elements }")
    end
end

"""
    add_node!(model, ID, x, y, z)

Adds a new node to the model.
"""
function add_node!(model::Model, ID::Int, 
    x::Real, y::Real, z::Real)
    # Check if node with the provided ID already exist in the model:
    haskey(model.nodes, ID) && throw(ArgumentError(styled"Node with ID {emphasis:$(ID)} already exist in the model."))

    # Create a new node and add it to the model:
    model.nodes[ID] = Node(ID, x, y, z)

    # Return the updated model:
    return model
end

"""
    del_node!(model, ID)

Deletes a node from the model.
"""
function del_node!(model::Model, ID::Int)
    # Check if node with the provided ID exist in the model:
    !haskey(model.nodes, ID) && throw(ArgumentError(styled"Node with ID {emphasis:$(ID)} does not exist in the model."))

    # Delete the node from the model:
    delete!(model.nodes, ID)

    # Return the updated model:
    return model
end


"""
    add_material!(model, ID, E, υ)

Adds a new material to the model.
"""
function add_material!(model::Model, ID::Int, 
    E::Real, ν::Real, ρ::Real)
    # Check if material with the provided ID already exist in the model:
    haskey(model.materials, ID) && throw(ArgumentError(styled"Material with ID {emphasis:$(ID)} already exist in the model."))

    # Create a new material and add it to the model:
    model.materials[ID] = Material(ID, E, ν, ρ)

    # Return the updated model:
    return model
end

"""
    del_material!(model, ID)

Deletes a material from the model.
"""
function del_material!(model::Model, ID::Int)
    # Check if material with the provided ID exist in the model:
    !haskey(model.materials, ID) && throw(ArgumentError(styled"Material with ID {emphasis:$(ID)} does not exist in the model."))

    # Delete the material from the model:
    delete!(model.materials, ID)

    # Return the updated model:
    return model
end

"""
    add_section!(model, ID, A, I_zz, I_yy, J)

Adds a new section to the model.
"""
function add_section!(model::Model, ID::Int, 
    A::Real, I_zz::Real, I_yy::Real, J::Real)
    # Check if section with the provided ID already exist in the model:
    haskey(model.sections, ID) && throw(ArgumentError(styled"Section with ID {emphasis:$(ID)} already exists in the model."))

    # Create a new section and add it to the model:
    model.sections[ID] = Section(ID, A, I_zz, I_yy, J)

    # Return the updated model:
    return model
end

"""
    del_section!(model, ID)

Deletes a section from the model.
"""
function del_section!(model::Model, ID::Int)
    # Check if section with the provided ID exist in the model:
    !haskey(model.sections, ID) && throw(ArgumentError(styled"Section with ID {emphasis:$(ID)} does not exist in the model."))

    # Delete the section from the model:
    delete!(model.sections, ID)

    # Return the updated model:
    return model
end

"""
    add_element!(model, ID, node_i_ID, node_j_ID, material_ID, section_ID; 
        ω = 0, 
        releases_i = [false, false, false, false, false, false],
        releases_j = [false, false, false, false, false, false])

Adds a new element to the model.
"""
function add_element!(model::Model, ID::Int,
    node_i_ID::Int, node_j_ID::Int, material_ID::Int, section_ID::Int;
    ω::Real = 0,
    releases_i::Vector{Bool} = [false, false, false, false, false, false],
    releases_j::Vector{Bool} = [false, false, false, false, false, false])
    # Check if element with the provided ID already exist in the model:
    haskey(model.elements, ID) && throw(ArgumentError(styled"Element with ID {emphasis:$(ID)} already exist in the model."))

    # Check if nodes with the provided IDs exist in the model:
    !haskey(model.nodes, node_i_ID) && throw(ArgumentError(styled"Node with ID {emphasis:$(node_i_ID)} does not exist in the model."))
    !haskey(model.nodes, node_j_ID) && throw(ArgumentError(styled"Node with ID {emphasis:$(node_j_ID)} does not exist in the model."))

    # Check if material with the provided ID exist in the model:
    !haskey(model.materials, material_ID) && throw(ArgumentError(styled"Material with ID {emphasis:$(material_ID)} does not exist in the model."))

    # Check if section with the provided ID exist in the model:
    !haskey(model.sections, section_ID) && throw(ArgumentError(styled"Section with ID {emphasis:$(section_ID)} does not exist in the model."))

    # Extract the information about the nodes of the element:
    node_i = model.nodes[node_i_ID]
    node_j = model.nodes[node_j_ID]
    x_i, y_i, z_i = node_i.x, node_i.y, node_i.z
    x_j, y_j, z_j = node_j.x, node_j.y, node_j.z

    # Extract the information about the material of the element:
    material = model.materials[material_ID]
    E, ν, ρ = material.E, material.ν, material.ρ

    # Extract the information about the section of the element:
    section = model.sections[section_ID]
    A, I_zz, I_yy, J = section.A, section.I_zz, section.I_yy, section.J

    # Compute the length of the element:
    L = _compute_L(x_i, y_i, z_i, x_j, y_j, z_j)

    # Compute the transformation matrices:
    γ, T = _compute_T(x_i, y_i, z_i, x_j, y_j, z_j, L, ω)

    # Compute the element elastic stiffness matrix in the local coordinate system:
    k_e_l = _compute_k_e_l(E, ν, A, I_zz, I_yy, J, L)

    # Compute the element elastic stiffness matrix in the global coordinate system:
    k_e_g = T' * k_e_l * T

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, k_e_g, k_e_g)

    # Compute the element geometric stiffness matrix in its local coordinate system (without the axial load term):
    k_g_l = _compute_k_g_l(A, I_zz, I_yy, L)

    # Compute the element geometric stiffness matrix in the global coordinate system (without the axial load term):
    k_g_g = T' * k_g_l * T

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, k_g_g, k_g_g)

    # Create a new element and add it to the model:
    model.elements[ID] = Element(ID,    
        node_i_ID, x_i, y_i, z_i, 
        node_j_ID, x_j, y_j, z_j, 
        E, ν, ρ, 
        A, I_zz, I_yy, J, 
        ω, releases_i, releases_j,
        L, γ, T, k_e_l, k_e_g, k_g_l, k_g_g)

    # Return the updated model:
    return model
end

"""
    del_element!(model, ID)

Deletes an element from the model.
"""
function del_element!(model::Model, ID::Int)
    # Check if element with the provided ID exist in the model:
    !haskey(model.elements, ID) && throw(ArgumentError(styled"Element with ID {emphasis:$(ID)} does not exist in the model."))

    # Delete the element from the model:
    delete!(model.elements, ID)

    # Return the updated model:
    return model
end

"""
    add_support!(model, ID, u_x, u_y, u_z, θ_x, θ_y, θ_z)

Adds a new support to the model.
"""
function add_support!(model::Model, ID::Int, 
    u_x::Bool, u_y::Bool, u_z::Bool,
    θ_x::Bool, θ_y::Bool, θ_z::Bool)
    # Check if the node with the provided ID exist in the model:
    !haskey(model.nodes, ID) && throw(ArgumentError(styled"Node with ID {emphasis:$(ID)} does not exist in the model."))

    # Check if the support at the node with the provided ID already exist in the model:
    haskey(model.supports, ID) && throw(ArgumentError(styled"Support at node with ID {emphasis:$(ID)} already exist in the model."))

    # Create a new support and add it to the model:
    model.supports[ID] = [u_x, u_y, u_z, θ_x, θ_y, θ_z]

    # Return the updated model:
    return model
end

"""
    del_support!(model, ID)

Deletes a support from the model.
"""
function del_support!(model::Model, ID::Int)
    # Check if the support at the node with the provided ID exist in the model:
    !haskey(model.supports, ID) && throw(ArgumentError(styled"Support at node with ID {emphasis:$(ID)} does not exist in the model."))

    # Delete the support from the model:
    delete!(model.supports, ID)

    # Return the updated model:
    return model
end

"""
    add_conc_load!(model, ID, F_x, F_y, F_z, M_x, M_y, M_z)

Adds a new concentrated load to the model.
"""
function add_conc_load!(model::Model, ID::Int,
    F_x::Real, F_y::Real, F_z::Real,
    M_x::Real, M_y::Real, M_z::Real)
    # Check if the node with the provided ID exist in the model:
    !haskey(model.nodes, ID) && throw(ArgumentError(styled"Node with ID {emphasis:$(ID)} does not exist in the model."))

    # Check if the concentrated load at the node with the provided ID already exist in the model:
    haskey(model.conc_loads, ID) && throw(ArgumentError(styled"Concentrated load at node with ID {emphasis:$(ID)} already exist in the model."))

    # Create a new concentrated load and add it to the model:
    model.conc_loads[ID] = [F_x, F_y, F_z, M_x, M_y, M_z]

    # Return the updated model:
    return model
end

"""
    del_conc_load!(model, ID)

Deletes a concentrated load from the model.
"""
function del_conc_load!(model::Model, ID::Int)
    # Check if the concentrated load at the node with the provided ID exist in the model:
    !haskey(model.conc_loads, ID) && throw(ArgumentError(styled"Concentrated load at node with ID {emphasis:$(ID)} does not exist in the model."))

    # Delete the concentrated load from the model:
    delete!(model.conc_loads, ID)

    # Return the updated model:
    return model
end

"""
    add_dist_load!(model, ID, q_x, q_y, q_z)

Adds a new distributed load to the model.
"""
function add_dist_load!(model::Model, ID::Int,
    q_x::Real, q_y::Real, q_z::Real)
    # Check if the element with the provided ID exist in the model:
    !haskey(model.elements, ID) && throw(ArgumentError(styled"Element with ID {emphasis:$(ID)} does not exist in the model."))

    # Check if the distributed load at the element with the provided ID already exist in the model:
    haskey(model.dist_loads, ID) && throw(ArgumentError(styled"Distributed load at element with ID {emphasis:$(ID)} already exist in the model."))

    # Extract the element information:
    L = model.elements[ID].L # Length of the element
    T = model.elements[ID].T # Transformation matrix of the element

    # Create a new distributed load and add it to the model:
    model.dist_loads[ID] = [q_x, q_y, q_z]

    # Compute the element fixed-end force vector in its local coordinate system:
    p_l = _compute_p_l(q_x, q_y, q_z, L)

    # Compute the element fixed-end force vector in the global coordinate system:
    p_g = T * p_l

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, p_g, p_g)

    # Add the element fixed-end force vectors to the model:
    model.p_l[ID] = p_l
    model.p_g[ID] = p_g

    # Return the updated model:
    return model
end

"""
    del_dist_load!(model, ID)

Deletes a distributed load from the model.
"""
function del_dist_load!(model::Model, ID::Int)
    # Check if the distributed load at the element with the provided ID exist in the model:
    !haskey(model.dist_loads, ID) && throw(ArgumentError(styled"Distributed load at element with ID {emphasis:$(ID)} does not exist in the model."))

    # Delete the distributed load from the model:
    delete!(model.dist_loads, ID)

    # Delete the element fixed-end force vector in its local coordinate system from the model:
    delete!(model.p_l, ID)

    # Delete the element fixed-end force vector in the global coordinate system from the model:
    delete!(model.p_g, ID)

    # Return the updated model:
    return model
end