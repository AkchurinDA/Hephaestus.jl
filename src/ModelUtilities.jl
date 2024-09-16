function add_node!(model::Model, tag::Int,
    x::Real, y::Real, z::Real)
    # Check if the node with the same tag already exists in the model:
    haskey(model.nodes, tag) && throw(ArgumentError("Node with this tag already exists in the model."))

    # Create a new node:
    node = Node(tag, x, y, z)

    # Add the node to the model:
    model.nodes[tag] = node

    # Return the updated model:
    return model
end

function del_node!(model::Model, tag::Int)
    # Check if the node with the same tag already exists in the model:
    haskey(model.nodes, tag) || throw(ArgumentError("Node with this tag does not exist in the model."))

    # Remove the node from the model:
    delete!(model.nodes, tag)

    # Return the updated model:
    return model
end

function add_material!(model::Model, tag::Int,
    E::Real, υ::Real, ρ::Real)
    # Check if the material with the same tag already exists in the model:
    haskey(model.materials, tag) && throw(ArgumentError("Material with this tag already exists in the model."))

    # Create a new material:
    material = Material(tag, E, υ, ρ)

    # Add the material to the model:
    model.materials[tag] = material

    # Return the updated model:
    return model
end

function del_material!(model::Model, tag::Int)
    # Check if the material with the same tag already exists in the model:
    haskey(model.materials, tag) || throw(ArgumentError("Material with this tag does not exist in the model."))

    # Remove the material from the model:
    delete!(model.materials, tag)

    # Return the updated model:
    return model
end

function add_section!(model::Model, tag::Int,
    A::Real, I_zz::Real, I_yy::Real, J::Real)
    # Check if the section with the same tag already exists in the model:
    haskey(model.sections, tag) && throw(ArgumentError("Section with this tag already exists in the model."))

    # Create a new section:
    section = Section(tag, A, I_zz, I_yy, J)

    # Add the section to the model:
    model.sections[tag] = section

    # Return the updated model:
    return model
end

function del_section!(model::Model, tag::Int)
    # Check if the section with the same tag already exists in the model:
    haskey(model.sections, tag) || throw(ArgumentError("Section with this tag does not exist in the model."))

    # Remove the section from the model:
    delete!(model.sections, tag)

    # Return the updated model:
    return model
end

function add_element!(model::Model, tag::Int,
    node_i_ID::Int, node_j_ID::Int, material_ID::Int, section_ID::Int, ω::Real)
    # Check if the element with the same tag already exists in the model:
    haskey(model.elements, tag) && throw(ArgumentError("Element with this tag already exists in the model."))

    # Check if the nodes, material, and section exist in the model:
    haskey(model.nodes    , node_i_ID  ) || throw(ArgumentError("Node (i) with this tag does not exist in the model."))
    haskey(model.nodes    , node_j_ID  ) || throw(ArgumentError("Node (j) with this tag does not exist in the model."))
    haskey(model.materials, material_ID) || throw(ArgumentError("Material with this tag does not exist in the model."))
    haskey(model.sections , section_ID ) || throw(ArgumentError("Section with this tag does not exist in the model." ))

    # Get the coordinates of the nodes:
    x_i, y_i, z_i = model.nodes[node_i_ID].x, model.nodes[node_i_ID].y, model.nodes[node_i_ID].z
    x_j, y_j, z_j = model.nodes[node_j_ID].x, model.nodes[node_j_ID].y, model.nodes[node_j_ID].z

    # Get the material properties:
    E = model.materials[material_ID].E
    ν = model.materials[material_ID].ν

    # Get the section properties:
    A    = model.sections[section_ID].A
    I_zz = model.sections[section_ID].I_zz
    I_yy = model.sections[section_ID].I_yy
    J    = model.sections[section_ID].J

    # Compute the length of the element:
    L = _compute_L(x_i, y_i, z_i, x_j, y_j, z_j)

    # Compute the local-to-global transformation matrix of the element:
    T = _compute_T(x_i, y_i, z_i, x_j, y_j, z_j, ω, L)

    # Compute the element elastic stiffness matrix in the LCS:
    k_e_l = _compute_k_e_l(E, ν, A, I_zz, I_yy, J, L)

    # Compute the element elastic stiffness matrix in the GCS:
    k_e_g = T' * k_e_l * T

    # Remove small values if any:
    map!(x -> abs(x) < eps() ? 0 : x, k_e_g, k_e_g)

    # Create a new element:
    element = Element(tag, node_i_ID, node_j_ID, material_ID, section_ID, ω, L, T, k_e_l, k_e_g)

    # Add the element to the model:
    model.elements[tag] = element

    # Return the updated model:
    return model
end

function del_element!(model::Model, tag::Int)
    # Check if the element with the same tag already exists in the model:
    haskey(model.elements, tag) || throw(ArgumentError("Element with this tag does not exist in the model."))

    # Remove the element from the model:
    delete!(model.elements, tag)

    # Return the updated model:
    return model
end

function add_support!(model::Model, tag::Int,
    u_x::Bool, u_y::Bool, u_z::Bool, 
    θ_x::Bool, θ_y::Bool, θ_z::Bool)
    # TODO

    # Create a new support:
    support = [u_x, u_y, u_z, θ_x, θ_y, θ_z]

    # Add the support to the model:
    model.supports[tag] = support

    # Return the updated model:
    return model
end

function del_support!(model::Model, tag::Int)
    # TODO

    # Remove the support from the model:
    delete!(model.supports, tag)

    # Return the updated model:
    return model
end

function add_nodal_load!(model::Model, tag::Int,
    F_x::Real, F_y::Real, F_z::Real, 
    M_x::Real, M_y::Real, M_z::Real)
    # TODO

    # Create a new nodal load:
    nodal_load = [F_x, F_y, F_z, M_x, M_y, M_z]

    # Add the nodal load to the model:
    model.nodal_loads[tag] = nodal_load

    # Return the updated model:
    return model
end

function del_nodal_load!(model::Model, tag::Int)
    # TODO

    # Remove the nodal load from the model:
    delete!(model.nodal_loads, tag)

    # Return the updated model:
    return model
end