function getpartitionindices(model::Model{D}) where {D}
    # Initialize the partition indices:
    partitionindices = fill(false, 6 * length(model.nodes))

    # Assemble the partition indices for a 2D model:
    if D == 2
        for (i, node) in enumerate(model.nodes)
            # Extract the free DOFs of the node:
            u_x, u_y = node.u_x, node.u_y
            θ_z      = node.θ_z

            # Assemble the partition indices for the free DOFs:
            @inbounds partitionindices[6 * i - 5] = !u_x
            @inbounds partitionindices[6 * i - 4] = !u_y
            @inbounds partitionindices[6 * i - 3] = false
            @inbounds partitionindices[6 * i - 2] = false
            @inbounds partitionindices[6 * i - 1] = false
            @inbounds partitionindices[6 * i    ] = !θ_z
        end
    end

    # Assemble the partition indices for a 3D model:
    if D == 3
        for (i, node) in enumerate(model.nodes)
            # Extract the free DOFs of the node:
            u_x, u_y, u_z = node.u_x, node.u_y, node.u_z
            θ_x, θ_y, θ_z = node.θ_x, node.θ_y, node.θ_z

            # Assemble the partition indices for the free DOFs:
            @inbounds partitionindices[6 * i - 5] = !u_x
            @inbounds partitionindices[6 * i - 4] = !u_y
            @inbounds partitionindices[6 * i - 3] = !u_z
            @inbounds partitionindices[6 * i - 2] = !θ_x
            @inbounds partitionindices[6 * i - 1] = !θ_y
            @inbounds partitionindices[6 * i    ] = !θ_z
        end
    end

    # Return the partition indices:
    return partitionindices
end

function assemble_K_e(model::Model)
    # Initialize the global elastic stiffness matrix:
    T   = promote_type([gettype(element) for element in model.elements]...)
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
        k_e_g = element.state.k_e_g

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

function assemble_K_g(model::Model)
    # Initialize the global geometric stiffness matrix:
    T   = promote_type([gettype(element) for element in model.elements]...)
    K_g = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global geometric stiffness matrix:
    assemble_K_g!(K_g, model)

    # Return the global geometric stiffness matrix:
    return K_g
end

function assemble_K_g!(K_g::AbstractMatrix{<:Real}, model::Model)
    # Assemble the global geometric stiffness matrix:
    for element in model.elements
        # Extact the element stiffness matrix in the global coordinate system:
        k_g_g = element.state.k_g_g

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

function assemble_K_m(model::Model)
    # Initialize the global material stiffness matrix:
    T   = promote_type([gettype(element) for element in model.elements]...)
    K_m = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global material stiffness matrix:
    assemble_K_m!(K_m, model)

    # Return the global elastic stiffness matrix:
    return K_e
end

function assemble_K_m!(K_m::AbstractMatrix{<:Real}, model::Model)

end

function assemble_M(model::Model)
    # Initialize the global mass matrix:
    T = promote_type([gettype(element) for element in model.elements]...)
    M = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global mass matrix:
    assemble_M!(M, model)

    # Return the global mass matrix:
    return M
end

function assemble_M!(M::AbstractMatrix{<:Real}, model::Model)
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

function assemble_F_c(model::Model)
    if isempty(model.concloads)
        # Initialize the global load vector due to concentrated loads:
        F_c = zeros(6 * length(model.nodes))
    else
        # Initialize the global load vector due to concentrated loads:
        T   = promote_type([gettype(concload) for concload in model.concloads]...)
        F_c = zeros(T, 6 * length(model.nodes))

        for concload in model.concloads
            # Extract the concentrated load vector in the global coordinate system:
            F_x, F_y, F_z = concload.F_x, concload.F_y, concload.F_z
            M_x, M_y, M_z = concload.M_x, concload.M_y, concload.M_z

            # Assemble the concentrated load vector into the global load vector:
            index = findfirst(x -> x.ID == concload.ID, model.nodes)
            @inbounds F_c[6 * index - 5] += F_x
            @inbounds F_c[6 * index - 4] += F_y
            @inbounds F_c[6 * index - 3] += F_z
            @inbounds F_c[6 * index - 2] += M_x
            @inbounds F_c[6 * index - 1] += M_y
            @inbounds F_c[6 * index    ] += M_z
        end
    end

    # Return the global load vector:
    return F_c
end

function assemble_F_d(model::Model)
    if isempty(model.distloads)
        # Initialize the global load vector due to distributed loads:
        F_d = zeros(6 * length(model.nodes))
    else
        # Initialize the global load vector due to distributed loads:
        T   = promote_type([gettype(distload) for distload in model.distloads]...)
        F_d = zeros(T, 6 * length(model.nodes))

        for distload in model.distloads
            # Extract the fixed-end force vector in the global coordinate system:
            p = distload.p

            # Find the element to which the distributed load is applied:
            element = model.elements[findfirst(x -> x.ID == distload.ID, model.elements)]

            # Assemble the fixed-end force vector into the global load vector:
            index_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
            index_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
            @inbounds F_d[(6 * index_i - 5):(6 * index_i)] += p[ 1:6]
            @inbounds F_d[(6 * index_j - 5):(6 * index_j)] += p[7:12]
        end
    end

    # Return the global load vector:
    return F_d
end
