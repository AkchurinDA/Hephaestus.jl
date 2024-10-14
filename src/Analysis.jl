abstract type AbstractAnalysisType  end
abstract type AbstractAnalysisCache end

"""
    solve(model::Model, analysis::AbstractAnalysisType)
    
Solve the model using the specified analysis type.
"""
function solve end

include("Analysis Types/LinearElasticAnalysis.jl")
include("Analysis Types/NonlinearElasticAnalysis.jl")
include("Analysis Types/ElasticBucklingAnalysis.jl")
include("Analysis Types/FreeVibrationAnalysis.jl")

function _assemble_K_e(model::Model, mapping::Dict{Int, Int})
    # Compute the number of nodes in the model:
    num_nodes = length(model.nodes)

    # Preallocate the global elastic stiffness matrix:
    T   = promote_type([_get_k_e_g_type(element) for element in values(model.elements)]...)
    K_e = zeros(T, 6 * num_nodes, 6 * num_nodes)

    # Assemble the global elastic stiffness matrix:
    for element in values(model.elements)
        # Extract the internal node IDs:
        internal_node_i_tag = mapping[element.node_i_tag]
        internal_node_j_tag = mapping[element.node_j_tag]

        # Set the base indices:
        base_index_i = 6 * (internal_node_i_tag - 1)
        base_index_j = 6 * (internal_node_j_tag - 1)

        # Set the ranges:
        range_i = (base_index_i + 1):(base_index_i + 6)
        range_j = (base_index_j + 1):(base_index_j + 6)

        # Add the element elastic stiffness matrix to the global elastic stiffness matrix:
        K_e[range_i, range_i] += element.k_e_g[1:6 , 1:6 ]
        K_e[range_i, range_j] += element.k_e_g[1:6 , 7:12]
        K_e[range_j, range_i] += element.k_e_g[7:12, 1:6 ]
        K_e[range_j, range_j] += element.k_e_g[7:12, 7:12]
    end

    # Return the global elastic stiffness matrix:
    return K_e
end

function _assemble_K_g(model::Model, mapping::Dict{Int, Int}, P::AbstractVector{<:Real})
    # Compute the number of nodes in the model:
    num_nodes = length(model.nodes)

    # Preallocate the global geometric stiffness matrix:
    T   = float(promote_type(eltype(P), [eltype(element.k_g_g) for element in values(model.elements)]...))
    K_g = zeros(T, 6 * num_nodes, 6 * num_nodes)

    # Assemble the global geometric stiffness matrix:
    for (i, element) in enumerate(values(model.elements))
        # Extract the internal node IDs:
        internal_node_i_tag = mapping[element.node_i_tag]
        internal_node_j_tag = mapping[element.node_j_tag]

        # Set the base indices:
        base_index_i = 6 * (internal_node_i_tag - 1)
        base_index_j = 6 * (internal_node_j_tag - 1)

        # Set the ranges:
        range_i = (base_index_i + 1):(base_index_i + 6)
        range_j = (base_index_j + 1):(base_index_j + 6)

        # Add the element geometric stiffness matrix to the global geometric stiffness matrix:
        K_g[range_i, range_i] += P[i] * element.k_g_g[1:6 , 1:6 ]
        K_g[range_i, range_j] += P[i] * element.k_g_g[1:6 , 7:12]
        K_g[range_j, range_i] += P[i] * element.k_g_g[7:12, 1:6 ]
        K_g[range_j, range_j] += P[i] * element.k_g_g[7:12, 7:12]
    end

    # Return the global geometric stiffness matrix:
    return K_g
end

function _assemble_M(model::Model, mapping::Dict{Int, Int})
    # Compute the number of nodes in the model:
    num_nodes = length(model.nodes)

    # Preallocate the global mass matrix:
    T = promote_type([_get_m_g_type(element) for element in values(model.elements)]...)
    M = zeros(T, 6 * num_nodes, 6 * num_nodes)

    # Assemble the global mass matrix:
    for element in values(model.elements)
        # Extract the internal node IDs:
        internal_node_i_tag = mapping[element.node_i_tag]
        internal_node_j_tag = mapping[element.node_j_tag]

        # Set the base indices:
        base_index_i = 6 * (internal_node_i_tag - 1)
        base_index_j = 6 * (internal_node_j_tag - 1)

        # Set the ranges:
        range_i = (base_index_i + 1):(base_index_i + 6)
        range_j = (base_index_j + 1):(base_index_j + 6)

        # Add the element mass matrix to the global mass matrix:
        M[range_i, range_i] += element.m_g[1:6 , 1:6 ]
        M[range_i, range_j] += element.m_g[1:6 , 7:12]
        M[range_j, range_i] += element.m_g[7:12, 1:6 ]
        M[range_j, range_j] += element.m_g[7:12, 7:12]
    end

    # Return the global mass matrix:
    return M
end

function _assemble_F(model::Model, mapping::Dict{Int, Int})
    # Compute the number of nodes in the model:
    num_nodes = length(model.nodes)
    
    # Assemble the global force vector:
    if isempty(model.conc_loads)
        F = zeros(6 * num_nodes)
    else
        # Preallocate the global force vector:
        T = promote_type([eltype(conc_load) for conc_load in values(model.conc_loads)]...)
        F = zeros(T, 6 * length(model.nodes))

        # Assemble the global force vector:
        for tag in keys(model.nodes)
            if haskey(model.conc_loads, tag)
                # Extract the internal node ID:
                internal_node_tag = mapping[tag]

                # Set the base index:
                base_index = 6 * (internal_node_tag - 1)

                # Set the range:
                range = (base_index + 1):(base_index + 6)

                # Add the concentrated load to the global force vector:
                F[range] += model.conc_loads[tag]
            end
        end
    end

    # Return the global force vector:
    return F
end


function _assemble_P(model::Model, mapping::Dict{Int, Int})
    # Compute the number of nodes in the model:
    num_nodes = length(model.nodes)

    if isempty(model.p_g)
        P = zeros(6 * num_nodes)
    else
        # Preallocate the global fixed-end force vector:
        T = promote_type([eltype(p_g) for p_g in values(model.p_g)]...)
        P = zeros(T, 6 * length(model.nodes))

        # Assemble the global fixed-end force vector:
        for (tag, element) in model.elements
            if haskey(model.dist_loads, tag)
                # Extract the internal node IDs of the nodes of an element:
                internal_node_i_tag = mapping[element.node_i_tag]
                internal_node_j_tag = mapping[element.node_j_tag]        

                # Set the base indices:
                base_index_i = 6 * (internal_node_i_tag - 1)
                base_index_j = 6 * (internal_node_j_tag - 1)

                # Set the ranges:
                range_i = (base_index_i + 1):(base_index_i + 6)
                range_j = (base_index_j + 1):(base_index_j + 6)

                # Add the fixed-end forces to the global fixed-end force vector:
                P[range_i] += model.p_g[tag][1:6 ]
                P[range_j] += model.p_g[tag][7:12]
            end
        end
    end

    # Return the global fixed-end force vector:
    return P
end

function _get_partition_indices(model::Model, mapping::Dict{Int, Int})
    # Compute the number of nodes in the model:
    num_nodes = length(model.nodes)

    # Initialize the partition indices:
    indices_f = fill(false, 6 * num_nodes)
    indices_s = fill(false, 6 * num_nodes)

    # Iterate over the nodes in the model:
    for tag in keys(model.nodes)
        # Get the internal node tag:
        internal_tag = mapping[tag]

        # Set the base index:
        base_index = 6 * (internal_tag - 1)

        # Check if the node's DOFs are free or supported:
        if haskey(model.supports, tag)
            # Extract the support conditions:
            support = model.supports[tag]

            # Iterate over the DOFs:
            for i in 1:6
                if support[i] 
                    # Supported DOF:
                    indices_s[base_index + i] = true
                else
                    # Free DOF:
                    indices_f[base_index + i] = true
                end
            end
        end
    end

    # Return the partitioning indices:
    return indices_f, indices_s
end

"""
    get_node_u_g(solution_cache::AbstractAnalysisCache, tag::Int)

Extracts the displacement vector of a node in the global coordinate system.
"""
function get_node_u_g(solution_cache::AbstractAnalysisCache, tag::Int)
    # Extract the displacement vector of the node in the global coordinate system:
    u_g = _get_node_u_g(tag, solution_cache.mapping, solution_cache.U)

    # Return the displacement vector of the node in the global coordinate system:
    return u_g
end

function _get_node_u_g(tag::Int, mapping::Dict{Int, Int}, U::Vector{<:Real})
    # Extract the internal node tag of the node:
    internal_node_tag = mapping[tag]

    # Set the base index:
    base_index = 6 * (internal_node_tag - 1)

    # Set the range:
    range = (base_index + 1):(base_index + 6)

    # Extract the displacement vector of the node in the global coordinate system:
    u_g = U[range]

    # Return the displacement vector of the node in the global coordinate system:
    return u_g
end

"""
    get_element_u_l(model::Model, solution_cache::AbstractAnalysisCache, tag::Int)

Extracts the element displacement vector in the local coordinate system.
"""
function get_element_u_l(model::Model, solution_cache::AbstractAnalysisCache, tag::Int)
    # Extract the element displacement vector in the local coordinate system:
    u_l = _get_element_u_l(model, tag, solution_cache.mapping, solution_cache.U)

    # Return the element displacement vector in the local coordinate system:
    return u_l
end

function _get_element_u_l(model::Model, tag::Int, mapping::Dict{Int, Int}, U::Vector{<:Real})
    # Extract the displacements at the nodes of the element in the global coordinate system:
    node_i_u_g = _get_node_u_g(model.elements[tag].node_i_tag, mapping, U)
    node_j_u_g = _get_node_u_g(model.elements[tag].node_j_tag, mapping, U)

    # Combine the displacements at the nodes of the element in the global coordinate system:
    u_g = [node_i_u_g; node_j_u_g]

    # Extract the transformation matrix of the element:
    Γ = model.elements[tag].Γ

    # Transform the displacements to the local coordinate system:
    u_l = Γ * u_g

    # Return the element displacement vector in the local coordinate system:
    return u_l
end

"""
    get_element_f_l(model::Model, solution_cache::AbstractAnalysisCache, tag::Int)

Extract the element force vector in the local coordinate system.
"""
function get_element_f_l(model::Model, solution_cache::AbstractAnalysisCache, tag::Int)
    # Extract the element force vector in the local coordinate system:
    f_l = _get_element_f_l(model, tag, solution_cache.mapping, solution_cache.U)

    # Return the element force vector in the local coordinate system:
    return f_l
end

function _get_element_f_l(model::Model, tag::Int, mapping::Dict{Int, Int}, U::Vector{<:Real})
    u_l = _get_element_u_l(model, tag, mapping, U)

    # Extract the element elastic stiffness matrix in the local coordinate system:
    k_e_l = model.elements[tag].k_e_l

    # Compute the element force vector in the local coordinate system:
    f_l = k_e_l * u_l

    # Return the element force vector in the local coordinate system:
    return f_l
end