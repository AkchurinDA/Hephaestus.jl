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

    # Partition the DOFs into free and supported:
    indices_f, indices_s, U_s = _partition_U(model)

    # Solve the system of equations:
    K_ff = model.K_e[indices_f, indices_f]
    K_fs = model.K_e[indices_f, indices_s]
    F_f  = model.F[indices_f]
    display(K_fs)
    display(U_s)
    U_f  = K_ff \ (F_f - K_fs * U_s)

    # Solve for the reactions:
    K_sf = model.K_e[indices_s, indices_f]
    K_ss = model.K_e[indices_s, indices_s]
    F_s = K_sf * U_f + K_ss * U_s

    # Unpartition the global displacement vector:
    U = _unpartition_U(model, indices_f, indices_s, U_f, U_s)

    return U
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

function _partition_U(model::Model)
    # Preallocate:
    indices_f = Int[] # Free DOFs
    indices_s = Int[] # Supported DOFs

    # Preallocate:
    T   = promote_type([typeof(node.u_x_enforced) for node in values(model.nodes)]...)
    U_s = T[] # Supported DOFs

    # Find the free and supported DOFs:
    for node in values(model.nodes)
        ID_i = node.ID_i

        if !node.u_x_supported && node.u_x_enforced == 0 # If the DOF is free and has no enforced displacement
            push!(indices_f, 6 * ID_i - 5)
        elseif node.u_x_enforced != 0 # If the DOF has an enforced displacement
            push!(indices_s, 6 * ID_i - 5)
            push!(U_s, node.u_x_enforced)
        else # If the DOF is supported
            push!(indices_s, 6 * ID_i - 5)
            push!(U_s, 0)
        end

        if !node.u_y_supported && node.u_y_enforced == 0 # If the DOF is free and has no enforced displacement
            push!(indices_f, 6 * ID_i - 4)
        elseif node.u_y_enforced != 0 # If the DOF has an enforced displacement
            push!(indices_s, 6 * ID_i - 4)
            push!(U_s, node.u_y_enforced)
        else # If the DOF is supported
            push!(indices_s, 6 * ID_i - 4)
            push!(U_s, 0)
        end

        if !node.u_z_supported && node.u_z_enforced == 0 # If the DOF is free and has no enforced displacement
            push!(indices_f, 6 * ID_i - 3)
        elseif node.u_z_enforced != 0 # If the DOF has an enforced displacement
            push!(indices_s, 6 * ID_i - 3)
            push!(U_s, node.u_z_enforced)
        else # If the DOF is supported
            push!(indices_s, 6 * ID_i - 3)
            push!(U_s, 0)
        end

        if !node.θ_x_supported && node.θ_x_enforced == 0 # If the DOF is free and has no enforced displacement
            push!(indices_f, 6 * ID_i - 2)
        elseif node.θ_x_enforced != 0 # If the DOF has an enforced displacement
            push!(indices_s, 6 * ID_i - 2)
            push!(U_s, node.θ_x_enforced)
        else # If the DOF is supported
            push!(indices_s, 6 * ID_i - 2)
            push!(U_s, 0)
        end

        if !node.θ_y_supported && node.θ_y_enforced == 0 # If the DOF is free and has no enforced displacement
            push!(indices_f, 6 * ID_i - 1)
        elseif node.θ_y_enforced != 0 # If the DOF has an enforced displacement
            push!(indices_s, 6 * ID_i - 1)
            push!(U_s, node.θ_y_enforced)
        else # If the DOF is supported
            push!(indices_s, 6 * ID_i - 1)
            push!(U_s, 0)
        end

        if !node.θ_z_supported && node.θ_z_enforced == 0 # If the DOF is free and has no enforced displacement
            push!(indices_f, 6 * ID_i)
        elseif node.θ_z_enforced != 0 # If the DOF has an enforced displacement
            push!(indices_s, 6 * ID_i)
            push!(U_s, node.θ_z_enforced)
        else # If the DOF is supported
            push!(indices_s, 6 * ID_i)
            push!(U_s, 0)
        end
    end

    return indices_f, indices_s, U_s
end

function _unpartition_U(model::Model, indices_f::Vector{<:Int}, indices_s::Vector{<:Int}, U_f::Vector{FDT}, U_s::Vector{SDT}) where {FDT<:Real, SDT<:Real}
    # Preallocate:
    T = promote_type(FDT, SDT)
    U = zeros(T, 6 * length(model.nodes))

    # Assign the free DOFs:
    for (i, index) in enumerate(indices_f)
        U[index] = U_f[i]
    end

    # Assign the supported DOFs:
    for (i, index) in enumerate(indices_s)
        U[index] = U_s[i]
    end

    # Return the unpartitioned global displacement vector:
    return U
end