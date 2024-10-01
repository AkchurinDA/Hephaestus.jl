function _get_partitioning_indices(model::Model, internal_node_IDs::Dict{Int, Int})
    # Preallocate:
    indices_f = Int[] # List of free DOFs
    indices_s = Int[] # List of supported DOFs

    # Get the partitioning indices:
    for node in values(model.nodes)
        # Get the internal node ID:
        internal_node_ID = internal_node_IDs[node.ID]

        # Set the base index:
        base_index = 6 * (internal_node_ID - 1)

        # Check if the node's DOFs are free or supported:
        if haskey(model.supports, node.ID)
            for i in 1:6
                model.supports[node.ID][i] ? push!(indices_s, base_index + i) : push!(indices_f, base_index + i)
            end
        else
            append!(indices_f, (base_index + 1):(base_index + 6))
        end
    end

    # Return the partitioning indices:
    return indices_f, indices_s
end

function _assemble_K_e(model::Model, internal_node_IDs::Dict{Int, Int})
    # Preallocate the global elastic stiffness matrix:
    T   = promote_type([eltype(element.k_e_g) for element in values(model.elements)]...)
    K_e = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global elastic stiffness matrix:
    for element in values(model.elements)
        # Extract the internal node IDs:
        internal_node_i_ID = internal_node_IDs[element.node_i_ID]
        internal_node_j_ID = internal_node_IDs[element.node_j_ID]

        # Set the base indices:
        base_index_i = 6 * (internal_node_i_ID - 1)
        base_index_j = 6 * (internal_node_j_ID - 1)
 
        # Set the ranges:
        range_i = (base_index_i + 1):(base_index_i + 6)
        range_j = (base_index_j + 1):(base_index_j + 6)

        # Add the element elastic stiffness matrix to the global elastic stiffness matrix:
        @inbounds K_e[range_i, range_i] += element.k_e_g[1:6 , 1:6 ]
        @inbounds K_e[range_i, range_j] += element.k_e_g[1:6 , 7:12]
        @inbounds K_e[range_j, range_i] += element.k_e_g[7:12, 1:6 ]
        @inbounds K_e[range_j, range_j] += element.k_e_g[7:12, 7:12]
    end

    # Return the global elastic stiffness matrix:
    return K_e
end

function _assemble_K_g(model::Model, internal_node_IDs::Dict{Int, Int}, P::Vector{<:Real})
    # Preallocate the global geometric stiffness matrix:
    T   = float(promote_type(eltype(P), [eltype(element.k_g_g) for element in values(model.elements)]...))
    K_g = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global geometric stiffness matrix:
    for (i, element) in enumerate(values(model.elements))
        # Extract the internal node IDs:
        internal_node_i_ID = internal_node_IDs[element.node_i_ID]
        internal_node_j_ID = internal_node_IDs[element.node_j_ID]

        # Set the base indices:
        base_index_i = 6 * (internal_node_i_ID - 1)
        base_index_j = 6 * (internal_node_j_ID - 1)

        # Set the ranges:
        range_i = (base_index_i + 1):(base_index_i + 6)
        range_j = (base_index_j + 1):(base_index_j + 6)

        # Add the element geometric stiffness matrix to the global geometric stiffness matrix:
        @inbounds K_g[range_i, range_i] += P[i] * element.k_g_g[1:6 , 1:6 ]
        @inbounds K_g[range_i, range_j] += P[i] * element.k_g_g[1:6 , 7:12]
        @inbounds K_g[range_j, range_i] += P[i] * element.k_g_g[7:12, 1:6 ]
        @inbounds K_g[range_j, range_j] += P[i] * element.k_g_g[7:12, 7:12]
    end

    # Return the global geometric stiffness matrix:
    return K_g
end

function _assemble_F(model::Model, internal_node_IDs::Dict{Int, Int})
    if isempty(model.conc_loads)
        F = zeros(6 * length(model.nodes))
    else
        # Preallocate the global force vector:
        T = promote_type([eltype(concentrated_load) for concentrated_load in values(model.conc_loads)]...)
        F = zeros(T, 6 * length(model.nodes))

        # Assemble the global force vector:
        for node in values(model.nodes)
            if haskey(model.conc_loads, node.ID)
                # Extract the internal node ID:
                internal_node_ID = internal_node_IDs[node.ID]

                # Set the base index:
                base_index = 6 * (internal_node_ID - 1)

                # Set the range:
                range = (base_index + 1):(base_index + 6)

                # Add the concentrated load to the global force vector:
                @inbounds F[range] += model.conc_loads[node.ID]
            end
        end
    end

    # Return the global force vector:
    return F
end

function _assemble_P(model::Model, internal_node_IDs::Dict{Int, Int})
    if isempty(model.p_g)
        P = zeros(6 * length(model.nodes))
    else
        # Preallocate the global fixed-end force vector:
        T = promote_type([eltype(p_g) for p_g in values(model.p_g)]...)
        P = zeros(T, 6 * length(model.nodes))

        # Assemble the global fixed-end force vector:
        for element in values(model.elements)
            if haskey(model.dist_loads, element.ID)
                # Extract the internal node IDs of the nodes of an element:
                node_i_ID_i = internal_node_IDs[element.node_i_ID]
                node_j_ID_i = internal_node_IDs[element.node_j_ID]

                # Set the base indices:
                base_index_i = 6 * (node_i_ID_i - 1)
                base_index_j = 6 * (node_j_ID_i - 1)

                # Set the ranges:
                range_i = (base_index_i + 1):(base_index_i + 6)
                range_j = (base_index_j + 1):(base_index_j + 6)

                # Add the fixed-end forces to the global fixed-end force vector:
                @inbounds P[range_i] += model.p_g[element.ID][1:6 ]
                @inbounds P[range_j] += model.p_g[element.ID][7:12]
            end
        end
    end

    # Return the global fixed-end force vector:
    return P
end