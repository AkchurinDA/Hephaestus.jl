# --------------------------------------------------
# ADDING
# --------------------------------------------------
"""
    add_node!(model::Model, ID::Int, 
        x::Real, y::Real, z::Real)

Adds a node to the model.
"""
function add_node!(model::Model, ID::Int, 
    x::Real, y::Real, z::Real)
    # Check if the node already exists:
    if haskey(model.nodes, ID)
        error("Node with ID $ID already exists in the model.")
    end

    # Add the node to the model:
    model.nodes[ID] = Node(ID, 
        x, y, z)

    # Return the updated model:
    return model
end

"""
    add_material!(model::Model, ID::Int, 
        E::Real, ν::Real, ρ::Real)

Adds a material to the model.
"""
function add_material!(model::Model, ID::Int, 
    E::Real, ν::Real, ρ::Real)
    # Check if the material already exists:
    if haskey(model.materials, ID)
        error("Material with ID $(ID) already exists in the model.")
    end

    # Add the material to the model:
    model.materials[ID] = Material(ID, 
        E, ν, ρ)

    # Return the updated model:
    return model
end

"""
    add_section!(model::Model, ID::Int, 
        A::Real, I_zz::Real, I_yy::Real, J::Real)

Adds a section to the model.
"""
function add_section!(model::Model, ID::Int, 
    A::Real, I_zz::Real, I_yy::Real, J::Real)
    # Check if the section already exists:
    if haskey(model.sections, ID)
        error("Section with ID $(ID) already exists in the model.")
    end

    # Add the section to the model:
    model.sections[ID] = Section(ID, 
        A, I_zz, I_yy, J)

    # Return the updated model:
    return model
end

"""
    add_element!(model::Model, ID::Int,
        node_i_ID::Int, node_j_ID::Int, material_ID::Int, section_ID::Int;
        releases_i::Vector{Bool} = [false, false, false, false, false, false],
        releases_j::Vector{Bool} = [false, false, false, false, false, false],
        ω::Real = 0)

Adds an element to the model.
"""
function add_element!(model::Model, ID::Int,
    node_i_ID::Int, node_j_ID::Int, material_ID::Int, section_ID::Int;
    releases_i::Vector{Bool} = [false, false, false, false, false, false],
    releases_j::Vector{Bool} = [false, false, false, false, false, false],
    ω::Real = 0)
    # Check if the element already exists:
    if haskey(model.elements, ID)
        error("Element with ID $(ID) already exists in the model.")
    end

    # Extract the information about the nodes:
    if !haskey(model.nodes, node_i_ID)
        error("Node with ID $(node_i_ID) does not exist in the model.")
    end

    if !haskey(model.nodes, node_j_ID)
        error("Node with ID $(node_j_ID) does not exist in the model.")
    end

    x_i, y_i, z_i = model.nodes[node_i_ID].x, model.nodes[node_i_ID].y, model.nodes[node_i_ID].z
    x_j, y_j, z_j = model.nodes[node_j_ID].x, model.nodes[node_j_ID].y, model.nodes[node_j_ID].z

    # Extract the information about the material:
    if !haskey(model.materials, material_ID)
        error("Material with ID $(material_ID) does not exist in the model.")
    end

    E, ν, ρ = model.materials[material_ID].E, model.materials[material_ID].ν, model.materials[material_ID].ρ

    # Extract the information about the section:
    if !haskey(model.sections, section_ID)
        error("Section with ID $(section_ID) does not exist in the model.")
    end

    A, I_zz, I_yy, J = model.sections[section_ID].A, model.sections[section_ID].I_zz, model.sections[section_ID].I_yy, model.sections[section_ID].J

    # Add the element to the model:
    model.elements[ID] = Element(ID, 
        node_i_ID, node_j_ID, material_ID, section_ID,
        releases_i, releases_j, ω,
        x_i, y_i, z_i, 
        x_j, y_j, z_j,
        E, ν, ρ,
        A, I_zz, I_yy, J)

    # Return the updated model:
    return model
end

"""
    add_support!(model::Model, ID::Int, 
        u_x::Bool, u_y::Bool, u_z::Bool, 
        θ_x::Bool, θ_y::Bool, θ_z::Bool)

Adds a support to the model.
"""
function add_support!(model::Model, ID::Int, # Node ID
    u_x::Bool, u_y::Bool, u_z::Bool, 
    θ_x::Bool, θ_y::Bool, θ_z::Bool)
    # Check if the support already exists:
    if haskey(model.supports, ID)
        error("Support with ID $(ID) already exists in the model.")
    end

    # Check if the node exists in the model:
    if !haskey(model.nodes, ID)
        error("Node with ID $(ID) does not exist in the model.")
    end

    # Add the support to the model:
    model.supports[ID] = [u_x, u_y, u_z, θ_x, θ_y, θ_z]

    # Return the updated model:
    return model
end

"""
    add_conc_load!(model::Model, ID::Int,
        F_x::Real, F_y::Real, F_z::Real,
        M_x::Real, M_y::Real, M_z::Real)

Adds a concentrated load to the model.
"""
function add_conc_load!(model::Model, ID::Int,
    F_x::Real, F_y::Real, F_z::Real,
    M_x::Real, M_y::Real, M_z::Real)
    # Check if the node exists in the model:
    if !haskey(model.nodes, ID)
        error("Node with ID $(ID) does not exist in the model.")
    end

    # Check if the concentrated load already exists:
    if haskey(model.conc_loads, ID)
        error("Concentrated load at node with ID $(ID) already exists in the model.")
    end

    # Add the concentrated load to the model:
    model.conc_loads[ID] = [F_x, F_y, F_z, M_x, M_y, M_z]

    # Return the updated model:
    return model
end

"""
    add_dist_load!(model::Model, ID::Int,
        q_x::Real, q_y::Real, q_z::Real)

Adds a distributed load to the model.
"""
function add_dist_load!(model::Model, ID::Int,
    q_x::Real, q_y::Real, q_z::Real)
    # Check if the element exists in the model:
    if !haskey(model.elements, ID)
        error("Element with ID $(ID) does not exist in the model.")
    end

    # Check if the distributed load already exists:
    if haskey(model.dist_loads, ID)
        error("Distributed load on element with ID $(ID) already exists in the model.")
    end

    # Extract the information about the element:
    L = model.elements[ID].L
    T = model.elements[ID].T

    # Add the distributed load to the model:
    model.dist_loads[ID] = [q_x, q_y, q_z]

    # Compute the fixed-end forces in the local coordinate system:
    p_l = _compute_p_l(q_x, q_y, q_z, L)

    # Transform the fixed-end forces to the global coordinate system:
    p_g = T * p_l

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, p_g, p_g)

    # Add the fixed-end forces to the model:
    model.p_l[ID] = p_l
    model.p_g[ID] = p_g

    # Return the updated model:
    return model
end

# --------------------------------------------------
# DELETING
# --------------------------------------------------
"""
    del_node!(model::Model, ID::Int)

Deletes a node from the model and all related elements, supports, and concentrated loads.
"""
function del_node!(model::Model, ID::Int)
    # Check if the node exists in the model:
    if !haskey(model.nodes, ID)
        error("Node with ID $(ID) does not exist in the model.")
    end

    # Check if the node is related to any element:
    related_elements = Int[]
    for element in values(model.elements)
        if element.node_i_ID == ID || element.node_j_ID == ID
            push!(related_elements, element.ID)
        end
    end

    if !isempty(related_elements)
        @warn "Node with ID $(ID) is related to elements with IDs $(related_elements). These elements will be deleted from the model."
        [del_element!(model, related_element) for related_element in related_elements]
    end

    # Check if the node is related to any support:
    related_supports = Int[]
    for key in keys(model.supports)
        if key == ID
            push!(related_supports, key)
        end
    end

    if !isempty(related_supports)
        @warn "Node with ID $(ID) is related to supports with IDs $(related_supports). These supports will be deleted from the model."
        [del_support!(model, related_support) for related_support in related_supports]
    end

    # Check if the node is related to any concentrated load:
    related_conc_loads = Int[]
    for key in keys(model.conc_loads)
        if key == ID
            push!(related_conc_loads, key)
        end
    end

    if !isempty(related_conc_loads)
        @warn "Node with ID $(ID) is related to concentrated loads with IDs $(related_conc_loads). These concentrated loads will be deleted from the model."
        [del_conc_load!(model, related_conc_load) for related_conc_load in related_conc_loads]
    end

    # Delete the node from the model:
    delete!(model.nodes, ID)

    # Return the updated model:
    return model
end

"""
    del_material!(model::Model, ID::Int)

Deletes a material from the model and all related elements.
"""
function del_material!(model::Model, ID::Int)
    # Check if the material exists in the model:
    if !haskey(model.materials, ID)
        error("Material with ID $(ID) does not exist in the model.")
    end

    # Check if the material is related to any element:
    related_elements = Int[]
    for element in values(model.elements)
        if element.material_ID == ID
            push!(related_elements, element.ID)
        end
    end

    if !isempty(related_elements)
        @warn "Material with ID $(ID) is related to elements with IDs $(related_elements). These elements will be deleted from the model."
        [del_element!(model, related_element) for related_element in related_elements]
    end

    # Delete the material from the model:
    delete!(model.materials, ID)

    # Return the updated model:
    return model
end

"""
    del_section!(model::Model, ID::Int)

Deletes a section from the model and all related elements.
"""
function del_section!(model::Model, ID::Int)
    # Check if the section exists in the model:
    if !haskey(model.sections, ID)
        error("Section with ID $(ID) does not exist in the model.")
    end

    # Check if the section is related to any element:
    related_elements = Int[]
    for element in values(model.elements)
        if element.section_ID == ID
            push!(related_elements, element.ID)
        end
    end

    if !isempty(related_elements)
        @warn "Section with ID $(ID) is related to elements with IDs $(related_elements). These elements will be deleted from the model."
        [del_element!(model, related_element) for related_element in related_elements]
    end

    # Delete the section from the model:
    delete!(model.sections, ID)

    # Return the updated model:
    return model
end

"""
    del_element!(model::Model, ID::Int)

Deletes an element from the model and all related distributed loads.
"""
function del_element!(model::Model, ID::Int)
    # Check if the element exists in the model:
    if !haskey(model.elements, ID)
        error("Element with ID $(ID) does not exist in the model.")
    end

    # Check if the element is related to any distributed load:
    related_dist_loads = Int[]
    for key in keys(model.dist_loads)
        if key == ID
            push!(related_dist_loads, key)
        end
    end

    if !isempty(related_dist_loads)
        @warn "Element with ID $(ID) is related to distributed loads with IDs $(related_dist_loads). These distributed loads will be deleted from the model."
        [del_dist_load!(model, related_dist_load) for related_dist_load in related_dist_loads]
        [delete!(model.p_l, related_dist_load) for related_dist_load in related_dist_loads]
        [delete!(model.p_g, related_dist_load) for related_dist_load in related_dist_loads]
    end

    # Delete the element from the model:
    delete!(model.elements, ID)

    # Return the updated model:
    return model
end

"""
    del_support!(model::Model, ID::Int)

Deletes a support from the model.
"""
function del_support!(model::Model, ID::Int)
    # Check if the support exists in the model:
    if !haskey(model.supports, ID)
        error("Support with ID $(ID) does not exist in the model.")
    end

    # Delete the support from the model:
    delete!(model.supports, ID)

    # Return the updated model:
    return model
end

"""
    del_conc_load!(model::Model, ID::Int)

Deletes a concentrated load from the model.
"""
function del_conc_load!(model::Model, ID::Int)
    # Check if the concentrated load exists in the model:
    if !haskey(model.conc_loads, ID)
        error("Concentrated load with ID $(ID) does not exist in the model.")
    end

    # Delete the concentrated load from the model:
    delete!(model.conc_loads, ID)

    # Return the updated model:
    return model
end

"""
    del_dist_load!(model::Model, ID::Int)

Deletes a distributed load from the model.
"""
function del_dist_load!(model::Model, ID::Int)
    # Check if the distributed load exists in the model:
    if !haskey(model.dist_loads, ID)
        error("Distributed load with ID $(ID) does not exist in the model.")
    end

    # Delete the distributed load from the model:
    delete!(model.dist_loads, ID)

    # Return the updated model:
    return model
end

# --------------------------------------------------
# EXTRACTING
# --------------------------------------------------
"""
    get_node_u_g(solution_cache::AbstractSolutionCache, ID::Int)

Extracts the displacement vector of a node in the global coordinate system.
"""
function get_node_u_g(solution_cache::AbstractSolutionCache, ID::Int)
    # Extract the internal node IDs:
    internal_node_IDs = solution_cache.internal_node_IDs

    # Extract the internal node ID of the node:
    internal_node_ID = internal_node_IDs[ID]

    # Set the base index:
    base_index = 6 * (internal_node_ID - 1)

    # Set the range:
    range = (base_index + 1):(base_index + 6)

    # Extract the displacement vector of the node in the global coordinate system:
    u_g = solution_cache.U[range]

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, u_g, u_g)

    # Return the displacement vector of the node in the global coordinate system:
    return u_g
end

"""
    get_element_u_l(model::Model, solution_cache::AbstractSolutionCache, ID::Int)

Extracts the element displacement vector in the local coordinate system.
"""
function get_element_u_l(model::Model, solution_cache::AbstractSolutionCache, ID::Int)
    # Extract the displacements at the nodes of the element in the global coordinate system:
    node_i_u_g = get_node_u_g(solution_cache, model.elements[ID].node_i_ID)
    node_j_u_g = get_node_u_g(solution_cache, model.elements[ID].node_j_ID)

    # Combine the displacements at the nodes of the element in the global coordinate system:
    u_g = [node_i_u_g; node_j_u_g]

    # Extract the transformation matrix of the element:
    Γ = model.elements[ID].Γ

    # Transform the displacements to the local coordinate system:
    u_l = Γ * u_g

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, u_l, u_l)

    # Return the element displacement vector in the local coordinate system:
    return u_l
end

"""
    get_element_f_l(model::Model, solution_cache::AbstractSolutionCache, ID::Int)

Extract the element force vector in the local coordinate system.
"""
function get_element_f_l(model::Model, solution_cache::AbstractSolutionCache, ID::Int)
    # Extract the element displacement vector in the local coordinate system:
    u_l = get_element_u_l(model, solution_cache, ID)

    # Extract the element elastic stiffness matrix in the local coordinate system:
    k_e_l = model.elements[ID].k_e_l

    # Compute the element force vector in the local coordinate system:
    f_l = k_e_l * u_l

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, f_l, f_l)

    # Return the element force vector in the local coordinate system:
    return f_l
end