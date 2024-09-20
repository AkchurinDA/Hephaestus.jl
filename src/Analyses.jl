abstract type AbstractAnalysis end

struct O1EAnalysis <: AbstractAnalysis

end

struct O1ECache
    K_e
    F
    P
    R
    U
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

    # Get the partition indices:
    indices_f, indices_s = _get_partition_indices(model, n_ID_mapping)

    # Assemble the global elastic stiffness matrix:
    K_e = _assemble_K_e(model, n_ID_mapping)
    K_e_ff = K_e[indices_f, indices_f]
    K_e_fs = K_e[indices_f, indices_s]
    K_e_sf = K_e[indices_s, indices_f]
    K_e_ss = K_e[indices_s, indices_s]

    # Assemble the global force vector:
    F = _assemble_F(model, n_ID_mapping)
    F_f = F[indices_f]
    F_s = F[indices_s]

    # Assemble the global fixed-end force vector:
    P = _assemble_P(model, n_ID_mapping)
    P_f = P[indices_f]
    P_s = P[indices_s]

    # Solve for the displacements at the unsupported DOFs:
    U_f = K_e_ff \ (F_f - P_f)

    # Assemble the global displacement vector:
    U = zeros(eltype(U_f), 6 * length(model.nodes))
    for (i, index) in enumerate(indices_f)
        U[index] = U_f[i]
    end

    # Compute the reaction forces at the supported DOFs:
    R_s = K_e_sf * U_f + P_s

    # Assemble the global reaction force vector:
    R = zeros(eltype(R_s), 6 * length(model.nodes))
    for (i, index) in enumerate(indices_s)
        R[index] = R_s[i]
    end

    # Return the displacements:
    return O1ECache(K_e, F, P, R, U)
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

function _assemble_F(model::Model, n_ID_mapping::Dict{Int, Int})
    if isempty(model.concentrated_loads)
        F = zeros(6 * length(model.nodes))
    else
        # Preallocate the global force vector:
        T = promote_type([eltype(concentrated_load) for concentrated_load in values(model.concentrated_loads)]...)
        F = zeros(T, 6 * length(model.nodes))

        for node in values(model.nodes)
            if haskey(model.concentrated_loads, node.ID)
                node_i_ID = n_ID_mapping[node.ID]

                range = (6 * (node_i_ID - 1) + 1):(6 * node_i_ID)

                @inbounds F[range] += model.concentrated_loads[node.ID]
            end
        end
    end

    # Return the global force vector:
    return F
end

function _assemble_P(model::Model, n_ID_mapping::Dict{Int, Int})
    if isempty(model.p_g)
        P = zeros(6 * length(model.nodes))
    else
        # Preallocate the global fixed-end force vector:
        T = promote_type([eltype(p_g) for p_g in values(model.p_g)]...)
        P = zeros(T, 6 * length(model.nodes))

        # Assemble the global fixed-end force vector:
        for element in values(model.elements)
            if haskey(model.distributed_loads, element.ID)
                # Extract the internal node IDs of the nodes of an element:
                node_i_ID_i = n_ID_mapping[element.node_i_ID]
                node_j_ID_i = n_ID_mapping[element.node_j_ID]

                range_i = (6 * (node_i_ID_i - 1) + 1):(6 * node_i_ID_i)
                range_j = (6 * (node_j_ID_i - 1) + 1):(6 * node_j_ID_i)

                @inbounds P[range_i] += model.p_g[element.ID][1:6 ]
                @inbounds P[range_j] += model.p_g[element.ID][7:12]
            end
        end
    end

    # Return the global fixed-end force vector:
    return P
end

function _get_partition_indices(model::Model, n_ID_mapping::Dict{Int, Int})
    # Preallocate:
    indices_f = Int[] # List of free DOFs
    indices_s = Int[] # List of supported DOFs

    # Get the partition indices:
    for node in values(model.nodes)
        node_tag_i = n_ID_mapping[node.ID]

        if haskey(model.supports, node.ID) # Check if the node's DOFs are fixed
            for i in 1:6
                model.supports[node.ID][i] ? push!(indices_s, 6 * (node_tag_i - 1) + i) : push!(indices_f, 6 * (node_tag_i - 1) + i)
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