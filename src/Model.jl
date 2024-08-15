include("Components/Node.jl")
include("Components/Material.jl")
include("Components/Section.jl")
include("Components/Element.jl")

@kwdef mutable struct Model
    "Vector that stores all the information about the nodes."
    nodes::Vector{Node} = Node[]
    "Vector that stores all the information about the materials."
    materials::Vector{Material} = Material[]
    "Vector that stores all the information about the sections."
    sections::Vector{Section} = Section[]
    "Vector that stores all the information about the elements."
    elements::Vector{Element} = Element[]

    "Vector that stores the user-defined IDs of the nodes."
    _node_IDs::Vector{Integer} = Integer[]
    "Vector that stores the user-defined IDs of the materials."
    _material_IDs::Vector{Integer} = Integer[]
    "Vector that stores the user-defined IDs of the sections."
    _section_IDs::Vector{Integer} = Integer[]
    "Vector that stores the user-defined IDs of the elements."
    _element_IDs::Vector{Integer} = Integer[]

    _num_nodes::Integer = 0
    _num_materials::Integer = 0
    _num_sections::Integer = 0
    _num_elements::Integer = 0
end

function addnode!(model::Model, ID::Integer, X::Real, Y::Real, Z::Real)
    # Check if the node ID is already in use:
    if ID in model._node_IDs
        error("Node ID $ID is already in use.")
    else
        # Increment the node counter:
        model._num_nodes += 1

        # Add the node ID to the list of node IDs:
        push!(model._node_IDs, ID)

        # Add the node to the model:
        push!(model.nodes, Node(ID, X, Y, Z))
    end

    # Return the model:
    return model
end

function addmaterial!(model::Model, ID::Integer, E::Real, v::Real)
    # Check if the material ID is already in use:
    if ID in model._material_IDs
        error("Material ID $ID is already in use.")
    else
        # Increment the material counter:
        model._num_materials += 1

        # Add the material ID to the list of material IDs:
        push!(model._material_IDs, ID)

        # Add the material to the model:
        push!(model.materials, Material(ID, E, v))
    end

    # Return the model:
    return model
end

function addsection!(model::Model, ID::Integer, A::Real, I_xx::Real, I_yy::Real, J::Real)
    # Check if the section ID is already in use:
    if ID in model._section_IDs
        error("Section ID $ID is already in use.")
    else
        # Increment the section counter:
        model._num_sections += 1
            
        # Add the section ID to the list of section IDs:
        push!(model._section_IDs, ID)

        # Add the section to the model:
        push!(model.sections, Section(ID, A, I_xx, I_yy, J))
    end

    # Return the model:
    return model
end

function addelement!(model::Model, ID::Integer, node_i_ID::Integer, node_j_ID::Integer, material_ID::Integer, section_ID::Integer)
    # Check if the element ID is already in use:
    if ID in model._element_IDs
        error("Element ID $ID is already in use.")
    else
        # Extract the nodes using the provided node IDs:
        node_i_index = findfirst(x -> x == node_i_ID, model._node_IDs)
        node_j_index = findfirst(x -> x == node_j_ID, model._node_IDs)

        # Check if the nodes exist:
        if isnothing(node_i_index) || isnothing(node_j_index)
            throw(ArgumentError("The element with ID $ID is referencing a non-existent node."))
        end

        # Extract the nodes:
        node_i = model.nodes[node_i_index]
        node_j = model.nodes[node_j_index]

        # Extract the material using the provided material ID:
        material_index = findfirst(x -> x == material_ID, model._material_IDs)

        # Check if the material exists:
        if isnothing(material_index)
            throw(ArgumentError("The element with ID $ID is referencing a non-existent material."))
        end

        # Extract the material:
        material = model.materials[material_index]

        # Extract the section using the provided section ID:
        section_index = findfirst(x -> x == section_ID, model._section_IDs)

        # Check if the section exists:
        if isnothing(section_index)
            throw(ArgumentError("The element with ID $ID is referencing a non-existent section."))
        end

        # Extract the section:
        section = model.sections[section_index]

        # Increment the element counter:
        model._num_elements += 1

        # Add the element ID to the list of element IDs:
        push!(model._element_IDs, ID)

        # Add the element to the model:
        push!(model.elements, Element(ID, node_i, node_j, material, section))
    end

    # Return the model:
    return model
end