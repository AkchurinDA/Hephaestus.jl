abstract type AbstractModel end

"""
    struct Model

A type representing the finite element model of a structure of interest.

$(FIELDS)
"""
@kwdef struct Model
    "Ordered dictionary that stores the nodes of the model"
    nodes         ::OrderedDict{Int, Node    } = OrderedDict{Int, Node    }()
    "Ordered dictionary that stores the materials of the model"
    materials     ::OrderedDict{Int, Material} = OrderedDict{Int, Material}()
    "Ordered dictionary that stores the materials of the model"
    sections      ::OrderedDict{Int, Section } = OrderedDict{Int, Section }()
    "Ordered dictionary that stores the elements of the model"
    elements      ::OrderedDict{Int, Element } = OrderedDict{Int, Element }()

    "Ordered dictionary that stores the supports of the model"
    supports      ::OrderedDict{Int, AbstractVector{Bool}} = OrderedDict{Int, AbstractVector{Bool}}()

    "Ordered dictionary that stores the concentrated loads of the model"
    conc_loads    ::OrderedDict{Int, AbstractVector{<:Real}} = OrderedDict{Int, AbstractVector{<:Real}}()
    "Ordered dictionary that stores the distributed loads of the model"
    dist_loads    ::OrderedDict{Int, AbstractVector{<:Real}} = OrderedDict{Int, AbstractVector{<:Real}}()

    "Ordered dictionary that stores the fixed-end forces in the local coordinate system"
    p_l           ::OrderedDict{Int, AbstractVector{<:Real}} = OrderedDict{Int, AbstractVector{<:Real}}()
    "Ordered dictionary that stores the fixed-end forces in the global coordinate system"
    p_g           ::OrderedDict{Int, AbstractVector{<:Real}} = OrderedDict{Int, AbstractVector{<:Real}}()
end

function Base.show(io::IO, model::Model)
    num_nodes      = length(model.nodes     )
    num_sections   = length(model.sections  )
    num_materials  = length(model.materials )
    num_elements   = length(model.elements  )
    num_supports   = length(model.supports  )
    num_conc_loads = length(model.conc_loads)
    num_dist_loads = length(model.dist_loads)

    if num_nodes == 0 && num_sections == 0 && num_materials == 0 && num_elements == 0 && num_supports == 0 && num_conc_loads == 0 && num_dist_loads == 0
        return println(io, styled"{cyan, bold: Empty model.}")
    else
        reserved_spacing = floor(Int, log10(maximum([num_nodes, num_materials, num_sections, num_elements, num_supports, num_conc_loads, num_dist_loads])))  + 1

        str_num_nodes      = lpad(num_nodes     , reserved_spacing)
        str_num_materials  = lpad(num_materials , reserved_spacing)
        str_num_sections   = lpad(num_sections  , reserved_spacing)
        str_num_elements   = lpad(num_elements  , reserved_spacing)
        str_num_supports   = lpad(num_supports  , reserved_spacing)
        str_num_conc_loads = lpad(num_conc_loads, reserved_spacing)
        str_num_dist_loads = lpad(num_dist_loads, reserved_spacing)

        println(io, styled"{cyan, bold: Model with:}")
        num_nodes      == 0 ? nothing : println(io, styled"{cyan: $(str_num_nodes     ) \t Nodes          }")
        num_sections   == 0 ? nothing : println(io, styled"{cyan: $(str_num_sections  ) \t Sections       }")
        num_materials  == 0 ? nothing : println(io, styled"{cyan: $(str_num_materials ) \t Materials      }")
        num_elements   == 0 ? nothing : println(io, styled"{cyan: $(str_num_elements  ) \t Elements       }")
        num_supports   == 0 ? nothing : println(io, styled"{cyan: $(str_num_supports  ) \t Supported nodes}")
        num_conc_loads == 0 ? nothing : println(io, styled"{cyan: $(str_num_conc_loads) \t Loaded nodes   }")
        num_dist_loads == 0 ? nothing : println(io, styled"{cyan: $(str_num_dist_loads) \t Loaded elements}")
    end

    return nothing
end

"""
    add_node!(model, tag, x, y, z)

Add a node to the model.
"""
function add_node!(model::Model, tag::Int,
    x::Real, y::Real, z::Real)
    # Check if the node already exists:
    if haskey(model.nodes, tag)
        throw(ArgumentError("Node with tag $(tag) already exists"))
    end

    # Add the node to the model:
    model.nodes[tag] = Node(x, y, z)

    # Return the updated model:
    return model
end

"""
    add_section!(model, tag, A, I_zz, I_yy, J)

Add a section to the model.
"""
function add_section!(model::Model, tag::Int,
    A::Real, I_zz::Real, I_yy::Real, J::Real)
    # Check if the section already exists:
    if haskey(model.sections, tag)
        throw(ArgumentError("Section with tag $(tag) already exists"))
    end

    # Add the section to the model:
    model.sections[tag] = Section(A, I_zz, I_yy, J)

    # Return the updated model:
    return model
end

"""
    add_material!(model, tag, E, ν, ρ)

Add a material to the model.
"""
function add_material!(model::Model, tag::Int,
    E::Real, ν::Real, ρ::Real)
    # Check if the material already exists:
    if haskey(model.materials, tag)
        throw(ArgumentError("Material with tag $(tag) already exists"))
    end

    # Add the material to the model:
    model.materials[tag] = Material(E, ν, ρ)

    # Return the updated model:
    return model
end

"""
    add_element!(model, tag, node_i_tag, node_j_tag, section_tag, material_tag)

Add an element to the model.    
"""
function add_element!(model::Model, tag::Int,
    node_i_tag::Int, node_j_tag::Int, section_tag::Int, material_tag::Int;
    ω::Real = 0, 
    releases_i::Vector{Bool} = [false, false, false, false, false, false],
    releases_j::Vector{Bool} = [false, false, false, false, false, false])
    # Check if the element already exists:
    if haskey(model.elements, tag)
        throw(ArgumentError("Element with tag $(tag) already exists"))
    end

    # Check if the node (i) exists:
    if !haskey(model.nodes, node_i_tag)
        throw(ArgumentError("Node with tag $(node_i_tag) does not exist"))
    end

    # Check if the node (j) exists:
    if !haskey(model.nodes, node_j_tag)
        throw(ArgumentError("Node with tag $(node_j_tag) does not exist"))
    end

    # Check if the section exists:
    if !haskey(model.sections, section_tag)
        throw(ArgumentError("Section with tag $(section_tag) does not exist"))
    end

    # Check if the material exists:
    if !haskey(model.materials, material_tag)
        throw(ArgumentError("Material with tag $(material_tag) does not exist"))
    end

    # Extract the nodal coordinates:
    node_i = model.nodes[node_i_tag]
    node_j = model.nodes[node_j_tag]
    x_i, y_i, z_i = node_i.x, node_i.y, node_i.z
    x_j, y_j, z_j = node_j.x, node_j.y, node_j.z

    # Extract the section properties:
    section = model.sections[section_tag]
    A, I_zz, I_yy, J = section.A, section.I_zz, section.I_yy, section.J

    # Extract the material properties:
    material = model.materials[material_tag]
    E, ν, ρ = material.E, material.ν, material.ρ

    # Add the element to the model:
    model.elements[tag] = Element(
        node_i_tag, x_i, y_i, z_i,
        node_j_tag, x_j, y_j, z_j,
        material_tag, E, ν, ρ,
        section_tag, A, I_zz, I_yy, J,
        ω,
        releases_i, releases_j)

    # Return the updated model:
    return model
end

"""
    add_support!(model, tag, u_x, u_y, u_z, θ_x, θ_y, θ_z)

Add a support to the model.
"""
function add_support!(model::Model, tag::Int, # Node tag
    u_x::Bool, u_y::Bool, u_z::Bool, 
    θ_x::Bool, θ_y::Bool, θ_z::Bool)
    # Check if the support already exists:
    if haskey(model.supports, tag)
        throw(ArgumentError("Support with tag $(tag) already exists"))
    end

    # Check if the node exists:
    if !haskey(model.nodes, tag)
        throw(ArgumentError("Node with tag $(tag) does not exist"))
    end

    # Add the support to the model:
    model.supports[tag] = [u_x, u_y, u_z, θ_x, θ_y, θ_z]

    # Return the updated model:
    return model
end

"""
    add_conc_load!(model, tag, F_x, F_y, F_z, M_x, M_y, M_z)

Add a concentrated load to the model.
"""
function add_conc_load!(model::Model, tag::Int, # Node tag
    F_x::Real, F_y::Real, F_z::Real, 
    M_x::Real, M_y::Real, M_z::Real)
    # Check if the concentrated load exists:
    if haskey(model.conc_loads, tag)
        throw(ArgumentError("Concentrated load with tag $(tag) already exists"))
    end

    # Check if the node exists:
    if !haskey(model.nodes, tag)
        throw(ArgumentError("Node with tag $(tag) does not exist"))
    end

    # Add the concentrated load to the model:
    model.conc_loads[tag] = [F_x, F_y, F_z, M_x, M_y, M_z]

    # Return the updated model:
    return model
end

"""
    add_dist_load!(model, tag, q_x, q_y, q_z)

Add a distributed load to the model.
"""
function add_dist_load!(model::Model, tag::Int, # Element tag
    q_x::Real, q_y::Real, q_z::Real)
    # Check if the distributed load exists:
    if haskey(model.dist_loads, tag)
        throw(ArgumentError("Distributed load with tag $(tag) already exists"))
    end

    # Check if the element exists:
    if !haskey(model.elements, tag)
        throw(ArgumentError("Element with tag $(tag) does not exist"))
    end

    # Add the distributed load to the model:
    model.dist_loads[tag] = [q_x, q_y, q_z]

    # Extract the information about the element:
    L = model.elements[tag].L
    Γ = model.elements[tag].Γ

    # Compute the fixed-end forces in the local coordinate system:
    p_l = _compute_p_l(q_x, q_y, q_z, L)

    # Transform the fixed-end forces to the global coordinate system:
    p_g = Γ * p_l

    # Add the fixed-end forces to the model:
    model.p_l[tag] = p_l
    model.p_g[tag] = p_g

    # Return the updated model:
    return model
end

