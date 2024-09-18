abstract type AbstractAnalysisType end

include("O1E.jl") # 1st-order elastic analysis
include("O2E.jl") # 2nd-order elastic analysis
include("EB.jl" ) # Elastic buckling analysis

function _get_partition_indices(model::Model, mapping::Dict{Int, Int})
    # Preallocate:
    indices_f = Int[] # List of free DOFs
    indices_s = Int[] # List of supported DOFs

    # Get the partition indices:
    for node in values(model.nodes)
        node_tag_i = mapping[node.tag]

        if haskey(model.supports, node.tag) # Check if the node's DOFs are fixed
            for i in 1:6
                model.supports[node.tag][i] ? push!(indices_s, 6 * (node_tag_i - 1) + i) : push!(indices_f, 6 * (node_tag_i - 1) + i)
            end
        else
            for i in 1:6
                push!(indices_f, 6 * (node_tag_i - 1) + i)
            end
        end
    end

    # Return the partition indices:
    return indices_f, indices_s
end

function _assemble_K_e(model::Model, mapping::Dict{Int, Int})
    # Preallocate the global elastic stiffness matrix:
    T   = promote_type([eltype(element.k_e_g) for element in values(model.elements)]...)
    K_e = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global elastic stiffness matrix:
    for element in values(model.elements)
        # Extract the internal node IDs of the nodes of an element:
        node_i_tag_i = mapping[element.node_i_tag]
        node_j_tag_i = mapping[element.node_j_tag]
 
        range_i = (6 * (node_i_tag_i - 1) + 1):(6 * node_i_tag_i)
        range_j = (6 * (node_j_tag_i - 1) + 1):(6 * node_j_tag_i)

        @inbounds K_e[range_i, range_i] += element.k_e_g[1:6 , 1:6 ]
        @inbounds K_e[range_i, range_j] += element.k_e_g[1:6 , 7:12]
        @inbounds K_e[range_j, range_i] += element.k_e_g[7:12, 1:6 ]
        @inbounds K_e[range_j, range_j] += element.k_e_g[7:12, 7:12]
    end

    # Return the global elastic stiffness matrix:
    return K_e
end

function _assemble_K_g(model::Model, mapping::Dict{Int, Int})
    # Preallocate the global elastic stiffness matrix:
    T   = promote_type([eltype(element.k_g_g) for element in values(model.elements)]...)
    K_g = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global elastic stiffness matrix:
    for element in values(model.elements)
        # Extract the internal node IDs of the nodes of an element:
        node_i_tag_i = mapping[element.node_i_tag]
        node_j_tag_i = mapping[element.node_j_tag]

        range_i = (6 * (node_i_tag_i - 1) + 1):(6 * node_i_tag_i)
        range_j = (6 * (node_j_tag_i - 1) + 1):(6 * node_j_tag_i)

        @inbounds K_g[range_i, range_i] += element.k_g_g[1:6 , 1:6 ]
        @inbounds K_g[range_i, range_j] += element.k_g_g[1:6 , 7:12]
        @inbounds K_g[range_j, range_i] += element.k_g_g[7:12, 1:6 ]
        @inbounds K_g[range_j, range_j] += element.k_g_g[7:12, 7:12]
    end

    # Return the global elastic stiffness matrix:
    return K_g
end

function _assemble_F(model::Model, mapping::Dict{Int, Int})
    # Preallocate the global force vector:
    T = promote_type([eltype(nodal_loads) for nodal_loads in values(model.nodal_loads)]...)
    F = zeros(T, 6 * length(model.nodes))

    # Assemble the global force vector:
    for node in values(model.nodes)
        if haskey(model.nodal_loads, node.tag)
            node_tag_i = mapping[node.tag]

            for i in 1:6
                if model.nodal_loads[node.tag][i] != 0
                    @inbounds F[6 * (node_tag_i - 1) + i] = model.nodal_loads[node.tag][i]
                end
            end
        end
    end

    # Return the global force vector:
    return F
end

function _unpartition_U(model::Model, indices_f::Vector{Int}, U_f::Vector{T}) where {T<:Real}
    # Preallocate the global displacement vector:
    U = zeros(T, 6 * length(model.nodes))

    # Assemble the global displacement vector:
    for (i, index) in enumerate(indices_f)
        U[index] = U_f[i]
    end

    # Return the global displacement vector:
    return U
end