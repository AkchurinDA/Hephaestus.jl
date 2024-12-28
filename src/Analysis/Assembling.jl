function getpartitionindices(model::Model)
    # Initialize the partition indices:
    partitionindices = fill(false, 6 * length(model.nodes))

    # Assemble the partition indices:
    if model.dimensionality == 2
        for (i, node) in enumerate(model.nodes)
            # Extract the free DOFs of the node:
            u_x, u_y = node.u_x, node.u_y
            θ_z      = node.θ_z

            # Assemble the partition indices for the free DOFs:
            partitionindices[6 * i - 5] = !u_x
            partitionindices[6 * i - 4] = !u_y
            partitionindices[6 * i - 3] = false
            partitionindices[6 * i - 2] = false
            partitionindices[6 * i - 1] = false
            partitionindices[6 * i    ] = !θ_z
        end
    elseif model.dimensionality == 3
        for (i, node) in enumerate(model.nodes)
            # Extract the free DOFs of the node:
            u_x, u_y, u_z = node.u_x, node.u_y, node.u_z
            θ_x, θ_y, θ_z = node.θ_x, node.θ_y, node.θ_z
    
            # Assemble the partition indices for the free DOFs:
            partitionindices[6 * i - 5] = !u_x
            partitionindices[6 * i - 4] = !u_y
            partitionindices[6 * i - 3] = !u_z
            partitionindices[6 * i - 2] = !θ_x
            partitionindices[6 * i - 1] = !θ_y
            partitionindices[6 * i    ] = !θ_z
        end
    end

    # Return the partition indices:
    return partitionindices
end

function assemble_K_e(model::Model)
    # Initialize the global elastic stiffness matrix:
    T   = promote_type([get_k_e_g_T(element) for element in model.elements]...)
    K_e = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global elastic stiffness matrix:
    assemble_K_e!(K_e, model)

    # Return the global elastic stiffness matrix:
    return K_e
end

function assemble_K_e!(K_e::AbstractMatrix{<:Real}, model::Model)
    # Assemble the global elastic stiffness matrix:
    for element in model.elements
        # Extact the element stiffness matrix in the global coordinate system:
        k_e_g = element.k_e_g

        # Assemble the element stiffness matrix into the global stiffness matrix:
        index_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        index_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds K_e[(6 * index_i - 5):(6 * index_i), (6 * index_i - 5):(6 * index_i)] += k_e_g[ 1:6,  1:6]
        @inbounds K_e[(6 * index_i - 5):(6 * index_i), (6 * index_j - 5):(6 * index_j)] += k_e_g[ 1:6, 7:12]
        @inbounds K_e[(6 * index_j - 5):(6 * index_j), (6 * index_i - 5):(6 * index_i)] += k_e_g[7:12,  1:6]
        @inbounds K_e[(6 * index_j - 5):(6 * index_j), (6 * index_j - 5):(6 * index_j)] += k_e_g[7:12, 7:12]
    end

    # Return the global elastic stiffness matrix:
    return K_e
end

function assemble_K_g(model::Model, N::AbstractVector{NT}) where {NT <: Real}
    # Make sure the length of the vector of axial loads is equal to the number of elements:
    @assert length(N) == length(model.elements) "The length of the vector of axial loads is not equal to the number of elements. Something is wrong."

    # Initialize the global geometric stiffness matrix:
    T   = promote_type([get_k_e_g_T(element) for element in model.elements]..., NT)
    K_g = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global geometric stiffness matrix:
    assemble_K_g!(K_g, model, N)

    # Return the global geometric stiffness matrix:
    return K_g
end

function assemble_K_g!(K_g::AbstractMatrix{<:Real}, model::Model, N::AbstractVector{<:Real})
    # Assemble the global geometric stiffness matrix:
    for (element, n) in zip(model.elements, N)
        # Extact the element stiffness matrix in the global coordinate system:
        k_g_g = n * element.k_g_g

        # Assemble the element stiffness matrix into the global stiffness matrix:
        index_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        index_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds K_g[(6 * index_i - 5):(6 * index_i), (6 * index_i - 5):(6 * index_i)] += k_g_g[ 1:6,  1:6]
        @inbounds K_g[(6 * index_i - 5):(6 * index_i), (6 * index_j - 5):(6 * index_j)] += k_g_g[ 1:6, 7:12]
        @inbounds K_g[(6 * index_j - 5):(6 * index_j), (6 * index_i - 5):(6 * index_i)] += k_g_g[7:12,  1:6]
        @inbounds K_g[(6 * index_j - 5):(6 * index_j), (6 * index_j - 5):(6 * index_j)] += k_g_g[7:12, 7:12]
    end

    # Return the global geometric stiffness matrix:
    return K_g
end

function assemble_M(model::Model)
    # Initialize the global mass matrix:
    T = promote_type([get_m_g_T(element) for element in model.elements]...)
    M = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global mass matrix:
    for element in model.elements
        # Extact the element mass matrix in the global coordinate system:
        m_g = element.m_g

        # Assemble the element mass matrix into the global mass matrix:
        index_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        index_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds M[(6 * index_i - 5):(6 * index_i), (6 * index_i - 5):(6 * index_i)] += m_g[ 1:6,  1:6]
        @inbounds M[(6 * index_i - 5):(6 * index_i), (6 * index_j - 5):(6 * index_j)] += m_g[ 1:6, 7:12]
        @inbounds M[(6 * index_j - 5):(6 * index_j), (6 * index_i - 5):(6 * index_i)] += m_g[7:12,  1:6]
        @inbounds M[(6 * index_j - 5):(6 * index_j), (6 * index_j - 5):(6 * index_j)] += m_g[7:12, 7:12]
    end

    # Return the global mass matrix:
    return M
end

function assemble_F_conc(model::Model)
    if isempty(model.concloads)
        # Initialize the global load vector due to concentrated loads:
        F_conc = zeros(6 * length(model.nodes))
    else
        # Initialize the global load vector due to concentrated loads:
        T = promote_type([get_concload_T(concload) for concload in model.concloads]...)
        F_conc = zeros(T, 6 * length(model.nodes))

        for concload in model.concloads
            # Extract the concentrated load vector in the global coordinate system:
            F_x, F_y, F_z = concload.F_x, concload.F_y, concload.F_z
            M_x, M_y, M_z = concload.M_x, concload.M_y, concload.M_z

            # Assemble the concentrated load vector into the global load vector:
            idx = findfirst(x -> x.ID == concload.ID, model.nodes)
            @inbounds F_conc[6 * idx - 5] += F_x
            @inbounds F_conc[6 * idx - 4] += F_y
            @inbounds F_conc[6 * idx - 3] += F_z
            @inbounds F_conc[6 * idx - 2] += M_x
            @inbounds F_conc[6 * idx - 1] += M_y
            @inbounds F_conc[6 * idx    ] += M_z
        end
    end

    # Return the global load vector:
    return F_conc
end

function assemble_F_dist(model::Model)
    if isempty(model.distloads)
        # Initialize the global load vector due to distributed loads:
        F_dist = zeros(6 * length(model.nodes))
    else
        # Initialize the global load vector due to distributed loads:
        T      = promote_type([get_distload_T(distload) for distload in model.distloads]...)
        F_dist = zeros(T, 6 * length(model.nodes))

        for distload in model.distloads
            # Extract the fixed-end force vector in the global coordinate system:
            p = distload.p

            # Find the element to which the distributed load is applied:
            element = model.elements[findfirst(x -> x.ID == distload.ID, model.elements)]

            # Assemble the fixed-end force vector into the global load vector:
            index_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
            index_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
            @inbounds F_dist[6 * index_i - 5] += p[ 1]
            @inbounds F_dist[6 * index_i - 4] += p[ 2]
            @inbounds F_dist[6 * index_i - 3] += p[ 3]
            @inbounds F_dist[6 * index_i - 2] += p[ 4]
            @inbounds F_dist[6 * index_i - 1] += p[ 5]
            @inbounds F_dist[6 * index_i    ] += p[ 6]
            @inbounds F_dist[6 * index_j - 5] += p[ 7]
            @inbounds F_dist[6 * index_j - 4] += p[ 8]
            @inbounds F_dist[6 * index_j - 3] += p[ 9]
            @inbounds F_dist[6 * index_j - 2] += p[10]
            @inbounds F_dist[6 * index_j - 1] += p[11]
            @inbounds F_dist[6 * index_j    ] += p[12]
        end
    end

    # Return the global load vector:
    return F_dist
end