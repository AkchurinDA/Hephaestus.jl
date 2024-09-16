function solve(model::Model, type::Symbol; log::Bool = true)
    # if log
    #     println(styled"{emphasis: -------------------------}")
    #     println(styled"{emphasis: Starting analaysis       }")
    #     println(styled"{emphasis: -------------------------}")
    # end

    # Reset the nodal displacements and rotations:
    for node in values(model.nodes)
        node.u_x = 0
        node.u_y = 0
        node.u_z = 0
        node.θ_x = 0
        node.θ_y = 0
        node.θ_z = 0
    end

    # Assign internal IDs to the nodes:
    for (i, node) in enumerate(values(model.nodes))
        node.ID_i = i
    end

    # Assign internal IDs to the elements:
    for (i, element) in enumerate(values(model.elements))
        element.ID_i = i
    end

    # Assemble the global stiffness matrix:
    model.K_e = _assemble_K_e(model)

    # Assemble the global force vector:
    model.F = _assemble_F(model)

    # Partition the global stiffness matrix:
    K_ff, K_fs, K_sf, K_ss = _partition_K(model, model.K_e)

    # Partition the global force vector:
    F_f, F_s = _partition_F(model, model.F)

    # Solve the system of equations:
    U_f = K_ff \ F_f
end

function _assemble_K_e(model::Model)
    # TODO: There should be a better way to do this.
    T   = promote_type([eltype(element.k_e_l_c) for element in values(model.elements)]...)
    K_e = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    for element in values(model.elements)
        # Extract the internal node IDs:
        node_i_ID_u = element.node_i_ID
        node_j_ID_u = element.node_j_ID
        node_i_ID_i = model.nodes[node_i_ID_u].ID_i
        node_j_ID_i = model.nodes[node_j_ID_u].ID_i

        @inbounds K_e[(6 * node_i_ID_i - 5):(6 * node_i_ID_i), (6 * node_i_ID_i - 5):(6 * node_i_ID_i)] += element.k_e_l_c[1:6 , 1:6 ]
        @inbounds K_e[(6 * node_i_ID_i - 5):(6 * node_i_ID_i), (6 * node_j_ID_i - 5):(6 * node_j_ID_i)] += element.k_e_l_c[1:6 , 7:12]
        @inbounds K_e[(6 * node_j_ID_i - 5):(6 * node_j_ID_i), (6 * node_i_ID_i - 5):(6 * node_i_ID_i)] += element.k_e_l_c[7:12, 1:6 ]
        @inbounds K_e[(6 * node_j_ID_i - 5):(6 * node_j_ID_i), (6 * node_j_ID_i - 5):(6 * node_j_ID_i)] += element.k_e_l_c[7:12, 7:12]
    end

    return K_e
end

function _assemble_F(model::Model)
    # TODO: There should be a better way to do this.
    T = promote_type([promote_type(
        typeof(node.F_x), typeof(node.F_y), typeof(node.F_z),
        typeof(node.M_x), typeof(node.M_y), typeof(node.M_z)) for node in values(model.nodes)]...)
    F = zeros(T, 6 * length(model.nodes))

    for node in values(model.nodes)
        @inbounds F[6 * node.ID_i - 5] += node.F_x
        @inbounds F[6 * node.ID_i - 4] += node.F_y 
        @inbounds F[6 * node.ID_i - 3] += node.F_z
        @inbounds F[6 * node.ID_i - 2] += node.M_x
        @inbounds F[6 * node.ID_i - 1] += node.M_y
        @inbounds F[6 * node.ID_i    ] += node.M_z
    end

    return F
end

function _partition_K(model::Model, K::Matrix{<:Real})
    # Preallocate:
    indices_f = Int[] # Free DOFs
    indices_s = Int[] # Supported DOFs

    # Find the free and supported DOFs:
    for node in values(model.nodes)
        ID_i = node.ID_i
        node.u_x_supported ? push!(indices_f, 6 * ID_i - 5) : push!(indices_s, 6 * ID_i - 5)
        node.u_y_supported ? push!(indices_f, 6 * ID_i - 4) : push!(indices_s, 6 * ID_i - 4)
        node.u_z_supported ? push!(indices_f, 6 * ID_i - 3) : push!(indices_s, 6 * ID_i - 3)
        node.θ_x_supported ? push!(indices_f, 6 * ID_i - 2) : push!(indices_s, 6 * ID_i - 2)
        node.θ_y_supported ? push!(indices_f, 6 * ID_i - 1) : push!(indices_s, 6 * ID_i - 1)
        node.θ_z_supported ? push!(indices_f, 6 * ID_i    ) : push!(indices_s, 6 * ID_i    )
    end

    # Partition the global stiffness matrix:
    K_ff = K[indices_f, indices_f]
    K_fs = K[indices_f, indices_s]
    K_sf = K[indices_s, indices_f]
    K_ss = K[indices_s, indices_s]

    # Return the partitioned stiffness matrix:
    return K_ff, K_fs, K_sf, K_ss
end

function _partition_F(model::Model, F::Vector{<:Real})
    # Preallocate:
    indices_f = Int[] # Free DOFs
    indices_s = Int[] # Supported DOFs

    # Find the free and supported DOFs:
    for node in values(model.nodes)
        ID_i = node.ID_i
        node.u_x_supported ? push!(indices_f, 6 * ID_i - 5) : push!(indices_s, 6 * ID_i - 5)
        node.u_y_supported ? push!(indices_f, 6 * ID_i - 4) : push!(indices_s, 6 * ID_i - 4)
        node.u_z_supported ? push!(indices_f, 6 * ID_i - 3) : push!(indices_s, 6 * ID_i - 3)
        node.θ_x_supported ? push!(indices_f, 6 * ID_i - 2) : push!(indices_s, 6 * ID_i - 2)
        node.θ_y_supported ? push!(indices_f, 6 * ID_i - 1) : push!(indices_s, 6 * ID_i - 1)
        node.θ_z_supported ? push!(indices_f, 6 * ID_i    ) : push!(indices_s, 6 * ID_i    )
    end

    # Partition the global force vector:
    F_f = F[indices_f]
    F_s = F[indices_s]

    # Return the partitioned force vector:
    return F_f, F_s
end