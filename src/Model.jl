"""
    struct Model

A type that represents the FE model of a structure.

To initialize a new model, use the [`Model()`](@ref) constructor.
"""
mutable struct Model
    "Dictionary that stores the information about the nodes of the model."
    nodes       ::OrderedDict{Int, Node    }
    "Dictionary that stores the information about the elements of the model."
    elements    ::OrderedDict{Int, Element }
    "Dictionary that stores the information about the materials of the model."
    materials   ::OrderedDict{Int, Material}
    "Dictionary that stores the information about the sections of the model."
    sections    ::OrderedDict{Int, Section }
    "Global elastic stiffness matrix of the model."
    K_e         ::Union{Nothing, Matrix}
    "Global force vector of the model."
    F           ::Union{Nothing, Vector}

    function Model()
        new(
            OrderedDict{Int, Node    }(),
            OrderedDict{Int, Element }(),
            OrderedDict{Int, Material}(),
            OrderedDict{Int, Section }(),
            nothing,
            nothing)
    end
end

# Pretty-printing the model:
function Base.show(io::IO, model::Model)
    if isempty(model.nodes) && isempty(model.materials) && isempty(model.sections) && isempty(model.elements)
        println(io, styled"{emphasis:A clean model with no nodes, materials, sections, or elements.}")
    else
        println(io, styled"{emphasis:Model of a structure with:}")
        println(io, styled"{emphasis:\t#\tComponent\tUser IDs}")
        !isempty(model.nodes    ) && println(io, styled"{grey:\t$(length(model.nodes    ))\tNodes    \t$(keys(model.nodes    ))}")
        !isempty(model.materials) && println(io, styled"{grey:\t$(length(model.materials))\tMaterials\t$(keys(model.materials))}")
        !isempty(model.sections ) && println(io, styled"{grey:\t$(length(model.sections ))\tSections \t$(keys(model.sections ))}")
        !isempty(model.elements ) && println(io, styled"{grey:\t$(length(model.elements ))\tElements \t$(keys(model.elements ))}")
    end
end

"""
    add_node!(model, ID, 
        x, y, z)

Adds a node to the model.
"""
function add_node!(model::Model, ID::Int,
    x::Real, y::Real, z::Real)
    # Check if the ID already exists:
    haskey(model.nodes, ID) && throw(ArgumentError("Node ID $ID already exists."))

    # Create the node:
    node = Node(ID, x, y, z)

    # Add the node to the model:
    model.nodes[ID] = node

    # Return the updated model:
    return model
end

"""
    del_node!(model, ID)

Deletes a node from the model.
"""
function del_node!(model::Model, ID::Int)
    # Check if the ID exists:
    !haskey(model.nodes, ID) && throw(ArgumentError("Node ID $ID does not exist."))

    # Delete the node from the model:
    delete!(model.nodes, ID)

    # Return the updated model:
    return model
end

"""
    add_material!(model, ID, 
        E, ν, ρ)

Adds a material to the model.
"""
function add_material!(model::Model, ID::Int,
    E::Real, ν::Real, ρ::Real)
    # Check if the ID already exists:
    haskey(model.materials, ID) && throw(ArgumentError("Material ID $ID already exists."))

    # Create the material:
    material = Material(ID, E, ν, ρ)

    # Add the material to the model:
    model.materials[ID] = material

    # Return the updated model:
    return model
end

"""
    del_material!(model, ID)

Deletes a material from the model.
"""
function del_material!(model::Model, ID::Int)
    # Check if the ID exists:
    !haskey(model.materials, ID) && throw(ArgumentError("Material ID $ID does not exist."))

    # Delete the material from the model:
    delete!(model.materials, ID)

    # Return the updated model:
    return model
end

"""
    add_section!(model, ID, 
        A, I_zz, I_yy, J)

Adds a section to the model.
"""
function add_section!(model::Model, ID::Int,
    A::Real, I_zz::Real, I_yy::Real, J::Real)
    # Check if the ID already exists:
    haskey(model.sections, ID) && throw(ArgumentError("Section ID $ID already exists."))

    # Create the section:
    section = Section(ID, A, I_zz, I_yy, J)

    # Add the section to the model:
    model.sections[ID] = section

    # Return the updated model:
    return model
end

"""
    del_section!(model, ID)

Deletes a section from the model.
"""
function del_section!(model::Model, ID::Int)
    # Check if the ID exists:
    !haskey(model.sections, ID) && throw(ArgumentError("Section ID $ID does not exist."))

    # Delete the section from the model:
    delete!(model.sections, ID)

    # Return the updated model:
    return model
end

"""
    add_element!(model, ID, 
        node_i_ID, node_j_ID, material_ID, section_ID; ω = 0)

Adds an element to the model.
"""
function add_element!(model::Model, ID::Int,
    node_i_ID::Int, node_j_ID::Int, material_ID::Int, section_ID::Int;
    ω::Real = 0, releases::Vector{Bool} = fill(false, 12))
    # Check if the ID already exists:
    haskey(model.elements, ID) && throw(ArgumentError("Element ID $ID already exists."))

    # Extract the nodes of the element:
    node_i = model.nodes[node_i_ID]
    node_j = model.nodes[node_j_ID]
    x_i, y_i, z_i = node_i.x, node_i.y, node_i.z
    x_j, y_j, z_j = node_j.x, node_j.y, node_j.z

    # Extract the material of the element:
    material = model.materials[material_ID]
    E, ν, ρ = material.E, material.ν, material.ρ

    # Extract the section of the element:
    section  = model.sections[section_ID]
    A, I_zz, I_yy, J = section.A, section.I_zz, section.I_yy, section.J

    # Compute the length of the element:
    L = _compute_L(x_i, y_i, z_i, x_j, y_j, z_j)

    # Compute the transformation matrix:
    T = _compute_T(x_i, y_i, z_i, x_j, y_j, z_j, L, ω)

    # Compute the uncondensed element elastic stiffness matrix in the local coordinate system:
    k_e_l = _compute_k_e_l(L, E, ν, A, I_zz, I_yy, J)

    # Compute the uncondensed element elastic stiffness matrix in the global coordinate system:
    k_e_g = T' * k_e_l * T

    # Compute the condensed element elastic stiffness matrix in the local coordinate system:
    k_e_l_c = _condense(k_e_l, releases)

    # Compute the condensed element elastic stiffness matrix in the global coordinate system:
    k_e_g_c = T' * k_e_l_c * T

    # Create the element:
    element = Element(ID, nothing,
        node_i_ID, x_i, y_i, z_i,
        node_j_ID, x_j, y_j, z_j,
        material_ID, E, ν, ρ,
        section_ID, A, I_zz, I_yy, J,
        ω, releases,
        L, T, k_e_l, k_e_g, k_e_l_c, k_e_g_c)

    # Add the element to the model:
    model.elements[ID] = element

    # Return the updated model:
    return model
end

"""
    del_element!(model, ID)

Deletes an element from the model.
"""
function del_element!(model::Model, ID::Int)
    # Check if the ID exists:
    !haskey(model.elements, ID) && throw(ArgumentError("Element ID $ID does not exist."))

    # Delete the element from the model:
    delete!(model.elements, ID)

    # Return the updated model:
    return model
end

"""
    add_support!(model, ID, 
        u_x, u_y, u_z, θ_x, θ_y, θ_z)

Fixes specified DOFs of a node in the model.
"""
function add_support!(model::Model, ID::Int,
    u_x_supported::Bool, u_y_supported::Bool, u_z_supported::Bool,
    θ_x_supported::Bool, θ_y_supported::Bool, θ_z_supported::Bool)
    # Check if the ID exists:
    !haskey(model.nodes, ID) && throw(ArgumentError("Node ID $ID does not exist."))

    # Check if the node already has supports applied to it:
    if model.nodes[ID].u_x_supported || model.nodes[ID].u_y_supported || model.nodes[ID].u_z_supported || 
       model.nodes[ID].θ_x_supported || model.nodes[ID].θ_y_supported || model.nodes[ID].θ_z_supported
        @warn "Node with ID = $ID already has supports applied to it. Overwriting them with the new values."
    end

    # Add the support to the node:
    model.nodes[ID].u_x_supported = u_x_supported
    model.nodes[ID].u_y_supported = u_y_supported
    model.nodes[ID].u_z_supported = u_z_supported
    model.nodes[ID].θ_x_supported = θ_x_supported
    model.nodes[ID].θ_y_supported = θ_y_supported
    model.nodes[ID].θ_z_supported = θ_z_supported

    # Return the updated model:
    return model
end

function del_support!(model::Model, ID::Int)
    # Check if the ID exists:
    !haskey(model.nodes, ID) && throw(ArgumentError("Node ID $ID does not exist."))

    # Remove the supports from the node:
    model.nodes[ID].u_x_supported = false
    model.nodes[ID].u_y_supported = false
    model.nodes[ID].u_z_supported = false
    model.nodes[ID].θ_x_supported = false
    model.nodes[ID].θ_y_supported = false
    model.nodes[ID].θ_z_supported = false

    # Return the updated model:
    return model
end

"""
    add_nodal_load!(model, ID, 
        F_x, F_y, F_z, 
        M_x, M_y, M_z)

Applies nodal loads to a node in the model.
"""
function add_nodal_load!(model::Model, ID::Int,
    F_x::NLT, F_y::NLT, F_z::NLT,
    M_x::NLT, M_y::NLT, M_z::NLT) where {NLT <: Real}
    # Check if the ID exists:
    !haskey(model.nodes, ID) && throw(ArgumentError("Node ID $ID does not exist."))

    # Check if the node already has nodal loads applied to it:
    if model.nodes[ID].F_x != 0 || model.nodes[ID].F_y != 0 || model.nodes[ID].F_z != 0 || 
       model.nodes[ID].M_x != 0 || model.nodes[ID].M_y != 0 || model.nodes[ID].M_z != 0
        @warn "Node with ID = $ID already has nodal loads applied to it. Overwriting them with the new values."
    end

    # Add the nodal load to the node:
    model.nodes[ID].F_x = F_x
    model.nodes[ID].F_y = F_y
    model.nodes[ID].F_z = F_z
    model.nodes[ID].M_x = M_x
    model.nodes[ID].M_y = M_y
    model.nodes[ID].M_z = M_z

    # Return the updated model:
    return model
end

add_nodal_load!(model::Model, ID::Int, 
    F_x::Real, F_y::Real, F_z::Real,
    M_x::Real, M_y::Real, M_z::Real) = 
    add_nodal_load!(model, ID, promote(F_x, F_y, F_z, M_x, M_y, M_z)...)

"""
    del_nodal_load!(model, ID)

Removes nodal loads from a node in the model.
"""
function del_nodal_load!(model::Model, ID::Int)
    # Check if the ID exists:
    !haskey(model.nodes, ID) && throw(ArgumentError("Node ID $ID does not exist."))

    # Set the nodal loads to zero:
    model.nodes[ID].F_x = 0
    model.nodes[ID].F_y = 0
    model.nodes[ID].F_z = 0
    model.nodes[ID].M_x = 0
    model.nodes[ID].M_y = 0
    model.nodes[ID].M_z = 0

    # Return the updated model:
    return model
end

function add_nodal_disp!(model::Model, ID::Int,
    u_x_enforced::Real, u_y_enforced::Real, u_z_enforced::Real,
    θ_x_enforced::Real, θ_y_enforced::Real, θ_z_enforced::Real)
    # Check if the ID exists:
    !haskey(model.nodes, ID) && throw(ArgumentError("Node ID $ID does not exist."))

    # Check if the node already has enforced displacements applied to it:
    if model.nodes[ID].u_x_enforced != 0 || model.nodes[ID].u_y_enforced != 0 || model.nodes[ID].u_z_enforced != 0 || 
       model.nodes[ID].θ_x_enforced != 0 || model.nodes[ID].θ_y_enforced != 0 || model.nodes[ID].θ_z_enforced != 0
        @warn "Node with ID = $ID already has enforced displacements applied to it. Overwriting them with the new values."
    end

    # Add the enforced displacements to the node:
    model.nodes[ID].u_x_enforced = u_x_enforced
    model.nodes[ID].u_y_enforced = u_y_enforced
    model.nodes[ID].u_z_enforced = u_z_enforced
    model.nodes[ID].θ_x_enforced = θ_x_enforced
    model.nodes[ID].θ_y_enforced = θ_y_enforced
    model.nodes[ID].θ_z_enforced = θ_z_enforced

    # Return the updated model:
    return model
end

add_nodal_disp!(model::Model, ID::Int, 
    u_x_enforced::Real, u_y_enforced::Real, u_z_enforced::Real,
    θ_x_enforced::Real, θ_y_enforced::Real, θ_z_enforced::Real) =
    add_nodal_disp!(model, ID, promote(u_x_enforced, u_y_enforced, u_z_enforced, θ_x_enforced, θ_y_enforced, θ_z_enforced)...)

function del_nodal_disp!(model::Model, ID::Int,
    u_x_enforced::NLT, u_y_enforced::NLT, u_z_enforced::NLT,
    θ_x_enforced::NLT, θ_y_enforced::NLT, θ_z_enforced::NLT) where {NLT <: Real}
    # Check if the ID exists:
    !haskey(model.nodes, ID) && throw(ArgumentError("Node ID $ID does not exist."))

    # Set the enforced displacements to zero:
    model.nodes[ID].u_x_enforced = 0
    model.nodes[ID].u_y_enforced = 0
    model.nodes[ID].u_z_enforced = 0
    model.nodes[ID].θ_x_enforced = 0
    model.nodes[ID].θ_y_enforced = 0
    model.nodes[ID].θ_z_enforced = 0

    # Return the updated model:
    return model
end

"""
    reset_model!(model)

Completely resets the model to a clean state.
"""
function reset_model!(model::Model)
    # Clear all dictionaries:
    empty!(model.nodes    )
    empty!(model.elements )
    empty!(model.materials)
    empty!(model.sections )

    # Return the updated model:
    return model
end