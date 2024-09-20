abstract type AbstractAnalysis end

struct O1EAnalysis <: AbstractAnalysis

end

struct O1ECache

end

struct O2EAnalysis <: AbstractAnalysis

end

struct O2ECache

end

function solve(model::Model, analysis::O1EAnalysis)
    # Assign internal IDs to nodes:
    n_ID_mapping = Dict{Int, Int}()
    for (i, node) in enumerate(values(model.nodes))
        n_ID_mapping[node.ID] = i
    end

    # Assign internal IDs to elements:
    e_ID_mapping = Dict{Int, Int}()
    for (i, element) in enumerate(values(model.elements))
        e_ID_mapping[element.ID] = i
    end
end

function _assemble_K_e(model::Model, n_ID_mapping::Dict{Int, Int})
    # Preallocate the global elastic stiffness matrix:
    T   = promote_type([eltype(element.k_e_g) for element in values(model.elements)]...)
    K_e = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global elastic stiffness matrix:
    for element in values(model.elements)
        # Extract the internal node IDs of the nodes of an element:
        node_i_ID_i = n_ID_mapping[element.node_i_ID]
        node_j_ID_i = n_ID_mapping[element.node_j_ID]
 
        range_i = (6 * (node_i_ID_i - 1) + 1):(6 * node_i_ID_i)
        range_j = (6 * (node_j_ID_i - 1) + 1):(6 * node_j_ID_i)

        @inbounds K_e[range_i, range_i] += element.k_e_g[1:6 , 1:6 ]
        @inbounds K_e[range_i, range_j] += element.k_e_g[1:6 , 7:12]
        @inbounds K_e[range_j, range_i] += element.k_e_g[7:12, 1:6 ]
        @inbounds K_e[range_j, range_j] += element.k_e_g[7:12, 7:12]
    end

    # Return the global elastic stiffness matrix:
    return K_e
end

function _assemble_K_g(model::Model, n_ID_mapping::Dict{Int, Int})
    # Preallocate the global geometric stiffness matrix:
    T   = promote_type([eltype(element.k_g_g) for element in values(model.elements)]...)
    K_g = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global geometric stiffness matrix:
    for element in values(model.elements)
        # Extract the internal node IDs of the nodes of an element:
        node_i_ID_i = n_ID_mapping[element.node_i_ID]
        node_j_ID_i = n_ID_mapping[element.node_j_ID]
 
        range_i = (6 * (node_i_ID_i - 1) + 1):(6 * node_i_ID_i)
        range_j = (6 * (node_j_ID_i - 1) + 1):(6 * node_j_ID_i)

        @inbounds K_g[range_i, range_i] += element.k_g_g[1:6 , 1:6 ]
        @inbounds K_g[range_i, range_j] += element.k_g_g[1:6 , 7:12]
        @inbounds K_g[range_j, range_i] += element.k_g_g[7:12, 1:6 ]
        @inbounds K_g[range_j, range_j] += element.k_g_g[7:12, 7:12]
    end

    # Return the global geometric stiffness matrix:
    return K_g
end