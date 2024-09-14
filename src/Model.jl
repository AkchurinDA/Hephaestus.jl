struct Node{T<:Real}
    "Unique identifier"
    ID          ::Int
    "``x``-coordinate"
    x           ::T
    "``y``-coordinate"
    y           ::T
    "``z``-coordinate"
    z           ::T
end

Node(ID::Int, x::Real, y::Real, z::Real) = Node(ID, promote(x, y, z)...)

struct Material{T<:Real}
    "Unique identifier"
    ID          ::Int
    "Young's modulus"
    E           ::T
    "Poisson's ratio"
    ν           ::T
end

Material(ID::Int, E::Real, ν::Real) = Material(ID, promote(E, ν)...)

struct Section{T<:Real}
    "Unique identifier"
    ID          ::Int
    "Area, ``A``"
    A           ::T
    "Moment of inertia about the local z-axis, ``I_{zz}``"
    I_zz        ::T
    "Moment of inertia about the local y-axis, ``I_{yy}``"
    I_yy        ::T
    "Polar moment of inertia, ``J``"
    J           ::T
end

Section(ID::Int, A::Real, I_zz::Real, I_yy::Real, J::Real) = Section(ID, promote(A, I_zz, I_yy, J)...)

@kwdef struct Element 
    "Unique identifier"
    ID          ::Int
    "Element's 1st node"
    node_1      ::Node
    "Element's 2nd node"
    node_2      ::Node
    "Element's material"
    material    ::Material
    "Element's section"
    section     ::Section
    "Element's orientation angle"
    ω           ::Real = 0
end

struct Model
    "List of model's nodes"
    nodes       ::Dict{Int, Node    }
    "List of model's materials"
    materials   ::Dict{Int, Material}
    "List of model's sections"
    sections    ::Dict{Int, Section }
    "List of model's elements"
    elements    ::Dict{Int, Element }

    # Initialize an empty model:
    function Model()
        new(
            Dict{Int, Node    }(), 
            Dict{Int, Material}(), 
            Dict{Int, Section }(), 
            Dict{Int, Element }())
    end
end

function add_node!(model::Model, ID::Int, x::Real, y::Real, z::Real)
    if haskey(model.nodes, ID)
        throw(ArgumentError("Node with ID = $ID already exists."))
    end

    model.nodes[ID] = Node(ID, x, y, z)

    return model
end

function delete_node!(model::Model, ID::Int)
    if !haskey(model.nodes, ID)
        throw(ArgumentError("Node with ID = $ID does not exist."))
    end

    delete!(model.nodes, ID)

    return model
end

function add_material!(model::Model, ID::Int, E::Real, ν::Real)
    if haskey(model.materials, ID)
        throw(ArgumentError("Material with ID = $ID already exists."))
    end

    model.materials[ID] = Material(ID, E, ν)

    return model
end

function delete_material!(model::Model, ID::Int)
    if !haskey(model.materials, ID)
        throw(ArgumentError("Material with ID = $ID does not exist."))
    end

    delete!(model.materials, ID)

    return model
end

function add_section!(model::Model, ID::Int, A::Real, I_zz::Real, I_yy::Real, J::Real)
    if haskey(model.sections, ID)
        throw(ArgumentError("Section with ID = $ID already exists."))
    end

    model.sections[ID] = Section(ID, A, I_zz, I_yy, J)

    return model
end

function delete_section!(model::Model, ID::Int)
    if !haskey(model.sections, ID)
        throw(ArgumentError("Section with ID = $ID does not exist."))
    end

    delete!(model.sections, ID)

    return model
end

function add_element!(model::Model, ID::Int, node_1::Node, node_2::Node, material::Material, section::Section)
    if haskey(model.elements, ID)
        throw(ArgumentError("Element with ID = $ID already exists."))
    end

    model.elements[ID] = Element(ID, node_1, node_2, material, section)

    return model
end

function delete_element!(model::Model, ID::Int)
    if !haskey(model.elements, ID)
        throw(ArgumentError("Element with ID = $ID does not exist."))
    end

    delete!(model.elements, ID)

    return model
end

function clear_model!(model::Model)
    # Empty the model's nodes, materials, sections, and elements:
    empty!(model.nodes    )
    empty!(model.materials)
    empty!(model.sections )
    empty!(model.elements )

    # Return the cleared model:
    return model
end

function compute_L(
    # Nodal coordinates of the element's 1st node:
    x_1::CT1, y_1::CT1, z_1::CT1,
    # Nodal coordinates of the element's 2nd node:
    x_2::CT2, y_2::CT2, z_2::CT2) where {CT1 <: Real, CT2 <: Real}
    # Compute the length of the element:
    L = sqrt((x_2 - x_1)^2 + (y_2 - y_1)^2 + (z_2 - z_1)^2)

    # Return the length of the element:
    return L
end

function compute_k_e(
    # Element properties:
    L::EPT,
    # Material properties:
    E::MPT, ν::MPT,
    # Section properties:
    A::SPT, I_zz::SPT, I_yy::SPT, J::SPT)::SparseMatrixCSC where {EPT <: Real, MPT <: Real, SPT <: Real}
    # Preallocate:
    T = promote_type(EPT, MPT, SPT)
    k_e = spzeros(T, 12, 12)

    # Compute the components of the element elastic stiffness matrix in its upper triangular part:
    @inbounds k_e[1 , 1 ] = +E * A / L
    @inbounds k_e[1 , 7 ] = -E * A / L
    @inbounds k_e[2 , 2 ] = +12 * E * I_zz / L^3
    @inbounds k_e[2 , 6 ] = +6 * E * I_zz / L^2
    @inbounds k_e[2 , 8 ] = -12 * E * I_zz / L^3
    @inbounds k_e[2 , 12] = +6 * E * I_zz / L^2
    @inbounds k_e[3 , 3 ] = +12 * E * I_yy / L^3
    @inbounds k_e[3 , 5 ] = -6 * E * I_yy / L^2
    @inbounds k_e[3 , 9 ] = -12 * E * I_yy / L^3
    @inbounds k_e[3 , 11] = -6 * E * I_yy / L^2
    @inbounds k_e[4 , 4 ] = +E * J / (2 * (1 + ν) * L)
    @inbounds k_e[4 , 10] = -E * J / (2 * (1 + ν) * L)
    @inbounds k_e[5 , 5 ] = +4 * E * I_yy / L
    @inbounds k_e[5 , 9 ] = +6 * E * I_yy / L^2
    @inbounds k_e[5 , 11] = +2 * E * I_yy / L
    @inbounds k_e[6 , 6 ] = +4 * E * I_zz / L
    @inbounds k_e[6 , 8 ] = -6 * E * I_zz / L^2
    @inbounds k_e[6 , 12] = +2 * E * I_zz / L
    @inbounds k_e[7 , 7 ] = +E * A / L
    @inbounds k_e[8 , 8 ] = -12 * E * I_zz / L^3
    @inbounds k_e[8 , 12] = -6 * E * I_zz / L^2
    @inbounds k_e[9 , 9 ] = +12 * E * I_yy / L^3
    @inbounds k_e[9 , 11] = +6 * E * I_yy / L^2
    @inbounds k_e[10, 10] = +E * J / (2 * (1 + ν) * L)
    @inbounds k_e[11, 11] = +12 * E * I_yy / L^3
    @inbounds k_e[12, 12] = +4 * E * I_zz / L

    # Compute the components of the element elastic stiffness matrix in its lower triangular part:
    for i in 1:12, j in (i + 1):12
        @inbounds k_e[j, i] = k_e[i, j]
    end

    # Drop zeros if any:
    dropzeros!(k_e)

    # Return the element elastic stiffness matrix:
    return k_e
end

function compute_k_g(
    # Element properties:
    L::EPT,
    # Section properties:
    A::SPT, I_zz::SPT, I_yy::SPT,
    # Axial load:
    P::ALT)::SparseMatrixCSC where {EPT <: Real, SPT <: Real, ALT <: Real}
    # Preallocate:
    T = promote_type(EPT, SPT, ALT)
    k_g = spzeros(T, 12, 12)

    # Compute the components of the element geometric stiffness matrix in its upper triangular part:
    @inbounds k_g[1 , 1 ] = +P / L
    @inbounds k_g[1 , 7 ] = -P / L
    @inbounds k_g[2 , 2 ] = +6 * P / (5 * L)
    @inbounds k_g[2 , 6 ] = +P / 10
    @inbounds k_g[2 , 8 ] = -6 * P / (5 * L)
    @inbounds k_g[2 , 12] = +P / 10
    @inbounds k_g[3 , 3 ] = +6 * P / (5 * L)
    @inbounds k_g[3 , 5 ] = -P / 10
    @inbounds k_g[3 , 9 ] = -6 * P / (5 * L)
    @inbounds k_g[3 , 11] = -P / 10
    @inbounds k_g[4 , 4 ] = +P * (I_zz + I_yy) / (A * L)
    @inbounds k_g[4 , 10] = -P * (I_zz + I_yy) / (A * L)
    @inbounds k_g[5 , 5 ] = +2 * P * L / 15
    @inbounds k_g[5 , 9 ] = +P / 10
    @inbounds k_g[5 , 11] = -P * L / 30
    @inbounds k_g[6 , 6 ] = +2 * P * L / 15
    @inbounds k_g[6 , 8 ] = -P / 10
    @inbounds k_g[6 , 12] = -P * L / 30
    @inbounds k_g[7 , 7 ] = +P / L
    @inbounds k_g[8 , 8 ] = +6 * P / (5 * L)
    @inbounds k_g[8 , 12] = -P / 10
    @inbounds k_g[9 , 9 ] = +6 * P / (5 * L)
    @inbounds k_g[9 , 11] = +P / 10
    @inbounds k_g[10, 10] = +P * (I_zz + I_yy) / (A * L)
    @inbounds k_g[11, 11] = +2 * P * L / 15
    @inbounds k_g[12, 12] = +2 * P * L / 15

    # Compute the components of the element geometric stiffness matrix in its lower triangular part:
    for i in 1:12, j in (i + 1):12
        @inbounds k_g[j, i] = k_g[i, j]
    end

    # Drop zeros if any:
    dropzeros!(k_g)

    # Return the element geometric stiffness matrix:
    return k_g
end

function assemble_K_e()

end