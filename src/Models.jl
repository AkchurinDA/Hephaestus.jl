@kwdef struct Model
    "Name of the model"
    name::String = "Model"
    "Nodes of the model"
    nodes::Vector{Node} = Vector{Node}()
    "Sections of the model"
    sections::Vector{Section} = Vector{Section}()
    "Materials of the model"
    materials::Vector{Material} = Vector{Material}()
    "Elements of the model"
    elements::Vector{Element} = Vector{Element}()
    "Concentrated loads of the model"
    concloads::Vector{ConcentratedLoad} = Vector{ConcentratedLoad}()
    "Distribution loads of the model"
    distloads::Vector{DistributedLoad} = Vector{DistributedLoad}()
end

function node!(model::Model, ID::Int, 
    x::Real, y::Real, z::Real; 
    u_x::Bool = false, u_y::Bool = false, u_z::Bool = false, 
    θ_x::Bool = false, θ_y::Bool = false, θ_z::Bool = false)
    # Check if the node already exists in the model:
    @assert ID ∉ getfield.(model.nodes, :ID) "Node already exists in the model."

    # Add the node to the model:
    push!(model.nodes, Node(ID, x, y, z, u_x, u_y, u_z, θ_x, θ_y, θ_z))

    # Return the updated model:
    return model
end

function section!(model::Model, ID::Int, A::Real, I_zz::Real, I_yy::Real, J::Real)
    # Check if the section already exists in the model:
    @assert ID ∉ getfield.(model.sections, :ID) "Section already exists in the model."

    # Add the section to the model:
    push!(model.sections, Section(ID, A, I_zz, I_yy, J))

    # Return the updated model:
    return model
end

function material!(model::Model, ID::Int, E::Real, ν::Real, ρ::Real)
    # Check if the material already exists in the model:
    @assert ID ∉ getfield.(model.materials, :ID) "Material already exists in the model."

    # Add the material to the model:
    push!(model.materials, Material(ID, E, ν, ρ))

    # Return the updated model:
    return model
end

function element!(model::Model, ID::Int, 
    node_i_ID::Int, node_j_ID::Int, 
    section_ID::Int, 
    material_ID::Int; 
    ω::Real = 0.0,
    releases_i::Vector{<:Bool} = [false, false, false, false, false, false],
    releases_j::Vector{<:Bool} = [false, false, false, false, false, false])
    # Check if the element already exists in the model:
    @assert ID ∉ getfield.(model.elements, :ID) "Element already exists in the model."

    # Check if the nodes exist in the model:
    @assert node_i_ID ∈ getfield.(model.nodes, :ID) "Node (``i``) with ID $(node_i_ID) does not exist in the model."
    @assert node_j_ID ∈ getfield.(model.nodes, :ID) "Node (``j``) with ID $(node_j_ID) does not exist in the model."

    # Check if the section exists in the model:
    @assert section_ID ∈ getfield.(model.sections, :ID) "Section with ID $(section_ID) does not exist in the model."

    # Check if the material exists in the model:
    @assert material_ID ∈ getfield.(model.materials, :ID) "Material with ID $(material_ID) does not exist in the model."

    # Extract the nodes:
    node_i = model.nodes[findfirst(x -> x.ID == node_i_ID, model.nodes)]
    node_j = model.nodes[findfirst(x -> x.ID == node_j_ID, model.nodes)]

    # Extract the section and material:
    section = model.sections[findfirst(x -> x.ID == section_ID, model.sections)]

    # Extract the material:
    material = model.materials[findfirst(x -> x.ID == material_ID, model.materials)]

    # Add the element to the model:
    push!(model.elements, Element(ID, node_i, node_j, section, material, ω, releases_i, releases_j))

    # Return the updated model:
    return model
end

function concload!(model::Model, ID::Int, 
    F_x::Real, F_y::Real, F_z::Real, 
    M_x::Real, M_y::Real, M_z::Real)
    # Check if the loads are zero:
    @assert F_x != 0 || F_y != 0 || F_z != 0 || M_x != 0 || M_y != 0 || M_z != 0 "All loads are zero. Aborting."

    # Check that the node exists in the model:
    @assert ID ∈ getfield.(model.nodes, :ID) "Node with ID $(ID) does not exist in the model."

    # Add the concentrated load to the model:
    push!(model.concloads, ConcentratedLoad(ID, F_x, F_y, F_z, M_x, M_y, M_z))

    # Return the updated model:
    return model
end

function distload!(model::Model, ID::Int, 
    w_x::Real, w_y::Real, w_z::Real; 
    cs::Symbol = :global)
    # Check if the loads are zero:
    @assert w_x != 0 || w_y != 0 || w_z != 0 "All loads are zero. Aborting."

    # Check that the element exists in the model:
    @assert ID ∈ getfield.(model.elements, :ID) "Element with ID $(element) does not exist in the model."

    # Check that a valid coordinate system is provided:
    @assert cs in [:local, :global] "Invalid coordinate system provided. Must be either `:local` or `:global`."

    # Add the distributed load to the model:
    push!(model.distloads, DistributedLoad(ID, w_x, w_y, w_z, cs))

    # Return the updated model:
    return model
end