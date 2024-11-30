@kwdef struct Model
    "Ordered dictionary that stores the nodes of the model"
    nodes    ::OrderedDict{Int, AbstractNode    } = OrderedDict{Int, AbstractNode    }()
    "Ordered dictionary that stores the materials of the model"
    materials::OrderedDict{Int, AbstractMaterial} = OrderedDict{Int, AbstractMaterial}()
    "Ordered dictionary that stores the sections of the model"
    sections ::OrderedDict{Int, AbstractSection } = OrderedDict{Int, AbstractSection }()
    "Ordered dictionary that stores the elements of the model"
    elements ::OrderedDict{Int, AbstractElement } = OrderedDict{Int, AbstractElement }()

    "Ordered dictionary that stores the supports applied to model"
    supports::OrderedDict{Int, Support} = OrderedDict{Int, Support}()

    "Ordered dictionary that stores the concentrated loads applied to the model"
    cloads::OrderedDict{Int, CLoad{<:Real}} = OrderedDict{Int, CLoad{<:Real}}()
    "Ordered dictionary that stores the distributed loads applied to the model"
    dloads::OrderedDict{Int, DLoad{<:Real}} = OrderedDict{Int, DLoad{<:Real}}()
end

# -------------------------
# ADDING
# -------------------------
function add_node!(model::Model, node_ID::Int, 
    x::Real, y::Real, z::Real)
    # Check if the node already exists in the model:
    !haskey(model.nodes, node_ID) || throw(ArgumentError("Node $(node_ID) already exists in the model."))

    # Add the node to the model:
    model.nodes[node_ID] = Node(x, y, z)

    # Return the updated model:
    return model
end

function add_material!(model::Model, material_ID::Int, 
    E::Real, ν::Real, ρ::Real)
    # Check if the material already exists in the model:
    !haskey(model.materials, material_ID) || throw(ArgumentError("Material $(material_ID) already exists in the model."))

    # Add the material to the model:
    model.materials[material_ID] = Material(E, ν, ρ)

    # Return the updated model:
    return model
end

function add_section!(model::Model, section_ID::Int, 
    A::Real, I_zz::Real, I_yy::Real, J::Real)
    # Check if the section already exists in the model:
    !haskey(model.sections, section_ID) || throw(ArgumentError("Section $(section_ID) already exists in the model."))

    # Add the section to the model:
    model.sections[section_ID] = Section(A, I_zz, I_yy, J)

    # Return the updated model:
    return model
end

function add_element!(model::Model, element_ID::Int, 
    node_i_ID::Int, 
    node_j_ID::Int,
    material_ID::Int, 
    section_ID::Int;
    ω::Real = 0)
    # Check if the element already exists in the model:
    !haskey(model.elements, element_ID) || throw(ArgumentError("Element $(element_ID) already exists in the model."))

    # Check if the nodes exist in the model:
    haskey(model.nodes, node_i_ID) || throw(ArgumentError("Node $(node_i_ID) does not exist in the model."))
    haskey(model.nodes, node_j_ID) || throw(ArgumentError("Node $(node_j_ID) does not exist in the model."))

    # Check if the material exists in the model:
    haskey(model.materials, material_ID) || throw(ArgumentError("Material $(material_ID) does not exist in the model."))

    # Check if the section exists in the model:
    haskey(model.sections, section_ID) || throw(ArgumentError("Section $(section_ID) does not exist in the model."))

    # Add the element to the model:
    model.elements[element_ID] = Element(
        # Node (i):
        node_i_ID, model.nodes[node_i_ID].x, model.nodes[node_i_ID].y, model.nodes[node_i_ID].z,
        # Node (j):
        node_j_ID, model.nodes[node_j_ID].x, model.nodes[node_j_ID].y, model.nodes[node_j_ID].z,
        # Material:
        material_ID, model.materials[material_ID].E, model.materials[material_ID].ν, model.materials[material_ID].ρ,
        # Section:
        section_ID, model.sections[section_ID].A, model.sections[section_ID].I_zz, model.sections[section_ID].I_yy, model.sections[section_ID].J,
        # Element orientation:
        ω)

    # Return the updated model:
    return model
end

function add_support!(model::Model, node_ID::Int,
    u_x::Bool, u_y::Bool, u_z::Bool, 
    θ_x::Bool, θ_y::Bool, θ_z::Bool)
    # Check if the support already exists:
    !haskey(model.supports, node_ID) || throw(ArgumentError("Support at node $(node_ID) already exists."))

    # Check if the node exists:
    haskey(model.nodes, node_ID) || throw(ArgumentError("Node $(node_ID) does not exist in the model."))

    # Add the support to the model:
    model.supports[node_ID] = Support(u_x, u_y, u_z, θ_x, θ_y, θ_z)

    # Return the updated model:
    return model
end

function add_cload!(model::Model, node_ID::Int,
    F_x::Real, F_y::Real, F_z::Real, 
    M_x::Real, M_y::Real, M_z::Real)
    # Check if the concentrated load exists:
    !haskey(model.cloads, node_ID) || throw(ArgumentError("Concentrated load with tag $(node_ID) already exists"))

    # Check if the node exists:
    haskey(model.nodes, node_ID) || throw(ArgumentError("Node with tag $(node_ID) does not exist"))

    # Add the concentrated load to the model:
    model.cloads[node_ID] = CLoad(F_x, F_y, F_z, M_x, M_y, M_z)

    # Return the updated model:
    return model
end

# -------------------------
# DELETING
# -------------------------
function del_node!(model::Model, node_ID::Int)
    # Check if the node exists in the model:
    haskey(model.nodes, node_ID) || throw(ArgumentError("Node $(node_ID) does not exist in the model."))

    # Check if the node is being used by any element:
    for (element_ID, element) in model.elements
        if (element.node_i_ID == node_ID || element.node_j_ID == node_ID)
            throw(ArgumentError("Cannot delete node $(node_ID) because it is being used by element $(element_ID)."))
        end
    end

    # Remove the node from the model:
    delete!(model.nodes, node_ID)

    # Return the updated model:
    return model
end

function del_material!(model::Model, mateiral_ID::Int)
    # Check if the material exists in the model:
    haskey(model.materials, mateiral_ID) || throw(ArgumentError("Material $(mateiral_ID) does not exist in the model."))

    # Remove the material from the model:
    delete!(model.materials, mateiral_ID)

    # Return the updated model:
    return model
end

function del_section!(model::Model, section_ID::Int)
    # Check if the section exists in the model:
    haskey(model.sections, section_ID) || throw(ArgumentError("Section $(section_ID) does not exist in the model."))

    # Remove the section from the model:
    delete!(model.sections, section_ID)

    # Return the updated model:
    return model
end

function del_element!(model::Model, element_ID::Int)
    # Check if the element exists in the model:
    haskey(model.elements, element_ID) || throw(ArgumentError("Element $(element_ID) does not exist in the model."))

    # Remove the element from the model:
    delete!(model.elements, element_ID)

    # Return the updated model:
    return model
end

function del_support!(model::Model, node_ID::Int)
    # Check if the support exists in the model:
    haskey(model.supports, node_ID) || throw(ArgumentError("Support at node $(node_ID) does not exist."))

    # Remove the support from the model:
    delete!(model.supports, node_ID)

    # Return the updated model:
    return model
end

function del_cload!(model::Model, node_ID::Int)
    # Check if the concentrated load exists in the model:
    haskey(model.conc_loads, node_ID) || throw(ArgumentError("Concentrated load at node $(node_ID) does not exist."))

    # Remove the concentrated load from the model:
    delete!(model.conc_loads, node_ID)

    # Return the updated model:
    return model
end