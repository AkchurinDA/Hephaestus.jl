function assemble_K_e(model::Model)
    # Initialize the global elastic stiffness matrix:
    K_e = zeros(Real, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global elastic stiffness matrix:
    for element in model.elements
        # Extact the transformation matrix of the element:
        Γ = element.Γ

        # Compute the element stiffness matrix in the global coordinate system:
        k_e_g = Γ' * element.k_e_l * Γ

        # Condense the element stiffness matrix if end releases are present:

        # Assemble the element stiffness matrix into the global stiffness matrix:
        idx_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        idx_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds K_e[(6 * idx_i - 5):(6 * idx_i), (6 * idx_i - 5):(6 * idx_i)] += k_e_g[ 1:6,  1:6]
        @inbounds K_e[(6 * idx_i - 5):(6 * idx_i), (6 * idx_j - 5):(6 * idx_j)] += k_e_g[ 1:6, 7:12]
        @inbounds K_e[(6 * idx_j - 5):(6 * idx_j), (6 * idx_i - 5):(6 * idx_i)] += k_e_g[7:12,  1:6]
        @inbounds K_e[(6 * idx_j - 5):(6 * idx_j), (6 * idx_j - 5):(6 * idx_j)] += k_e_g[7:12, 7:12]
    end

    # Return the global elastic stiffness matrix:
    return K_e
end

function assemble_K_g(model::Model, P::Vector{<:Real})
    # Make sure the length of the vector of axial loads is equal to the number of elements:
    @assert length(P) == length(model.elements) "The length of the vector of axial loads is not equal to the number of elements. Something is wrong."

    # Initialize the global geometric stiffness matrix:
    K_g = zeros(6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global geometric stiffness matrix:
    for (element, p) in zip(model.elements, P)
        # Extact the transformation matrix of the element:
        Γ = element.Γ

        # Compute the element stiffness matrix in the global coordinate system:
        k_e_g = p * Γ' * element.k_e_l * Γ

        # Condense the element stiffness matrix if end releases are present:

        # Assemble the element stiffness matrix into the global stiffness matrix:
        idx_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        idx_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds K_g[(6 * idx_i - 5):(6 * idx_i), (6 * idx_i - 5):(6 * idx_i)] += k_e_g[ 1:6,  1:6]
        @inbounds K_g[(6 * idx_i - 5):(6 * idx_i), (6 * idx_j - 5):(6 * idx_j)] += k_e_g[ 1:6, 7:12]
        @inbounds K_g[(6 * idx_j - 5):(6 * idx_j), (6 * idx_i - 5):(6 * idx_i)] += k_e_g[7:12,  1:6]
        @inbounds K_g[(6 * idx_j - 5):(6 * idx_j), (6 * idx_j - 5):(6 * idx_j)] += k_e_g[7:12, 7:12]
    end

    # Return the global geometric stiffness matrix:
    return K_g
end

function assemble_M(model::Model)
    # Initialize the global mass matrix:
    M = zeros(6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global mass matrix:
    for element in model.elements
        # Extact the transformation matrix of the element:
        Γ = element.Γ

        # Compute the element mass matrix in the global coordinate system:
        m_g = Γ' * element.m_l * Γ

        # Condense the element mass matrix if end releases are present:

        # Assemble the element mass matrix into the global mass matrix:
        idx_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        idx_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds M[(6 * idx_i - 5):(6 * idx_i), (6 * idx_i - 5):(6 * idx_i)] += m_g[ 1:6,  1:6]
        @inbounds M[(6 * idx_i - 5):(6 * idx_i), (6 * idx_j - 5):(6 * idx_j)] += m_g[ 1:6, 7:12]
        @inbounds M[(6 * idx_j - 5):(6 * idx_j), (6 * idx_i - 5):(6 * idx_i)] += m_g[7:12,  1:6]
        @inbounds M[(6 * idx_j - 5):(6 * idx_j), (6 * idx_j - 5):(6 * idx_j)] += m_g[7:12, 7:12]
    end
end

function assemble_F_conc(model::Model)
    # Initialize the global load vector due to concentrated loads:
    F_conc = zeros(Real, 6 * length(model.nodes))

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

    # Return the global load vector:
    return F_conc
end

function assemble_F_dist(model::Model)
    # Initialize the global load vector due to distributed loads:
    F_dist = zeros(Real, 6 * length(model.nodes))

    for distload in model.distloads
        # Find the element to which the distributed load is applied:
        element = model.elements[findfirst(x -> x.ID == distload.element.ID, model.elements)]

        # Extact the transformation matrix of the element:
        Γ = element.Γ

        # Extract the distributed load vector in the global coordinate system:
        w_x, w_y, w_z = distload.w_x, distload.w_y, distload.w_z
        
        # Extract the coordinate system in which the distributed loads are defined:
        cs = distload.cs

        # If the distributed loads are defined in the global coordinate system, convert them to the local coordinate system:
        if cs == :global 
            @error "Distributed loads defined in the global coordinate system are not yet supported."
            
            # Find the nodes of the element:
            node_i = model.nodes[findfirst(x -> x.ID == element.node_i.ID, model.nodes)]
            node_j = model.nodes[findfirst(x -> x.ID == element.node_j.ID, model.nodes)]

            # Extract the coordinates of the nodes:
            x_i, y_i, z_i = node_i.x, node_i.y, node_i.z
            x_j, y_j, z_j = node_j.x, node_j.y, node_j.z

            # Compute the element length projections:
            L_x = abs(x_j - x_i)
            L_y = abs(y_j - y_i)
            L_z = abs(z_j - z_i)

            # Extract the element length:
            L = element.L

            # Find the resultants of the distributed loads:
            R_x = w_x * L_x
            R_y = w_y * L_y
            R_z = w_z * L_z

            # TODO: Finish...
        end

        # Compute the fixed-end force vector:
        p_l = [
            -w_x * L / 2     ; # F_x_i
            -w_y * L / 2     ; # F_y_i
            -w_z * L / 2     ; # F_z_i
            0                ; # M_x_i
            -w_z * L ^ 2 / 12; # M_y_i
            -w_y * L ^ 2 / 12; # M_z_i
            +w_x * L / 2     ; # F_x_j
            -w_y * L / 2     ; # F_y_j
            -w_z * L / 2     ; # F_z_j
            0                ; # M_x_j
            +w_z * L ^ 2 / 12; # M_y_j
            +w_y * L ^ 2 / 12] # M_z_j

        # Convert the fixed-end force vector to the global coordinate system:
        p_g = Γ * p_l

        # Assemble the fixed-end force vector into the global load vector:
        idx_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
        idx_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)
        @inbounds F_dist[6 * idx_i - 5] += p_g[ 1]
        @inbounds F_dist[6 * idx_i - 4] += p_g[ 2]
        @inbounds F_dist[6 * idx_i - 3] += p_g[ 3]
        @inbounds F_dist[6 * idx_i - 2] += p_g[ 4]
        @inbounds F_dist[6 * idx_i - 1] += p_g[ 5]
        @inbounds F_dist[6 * idx_i    ] += p_g[ 6]
        @inbounds F_dist[6 * idx_j - 5] += p_g[ 7]
        @inbounds F_dist[6 * idx_j - 4] += p_g[ 8]
        @inbounds F_dist[6 * idx_j - 3] += p_g[ 9]
        @inbounds F_dist[6 * idx_j - 2] += p_g[10]
        @inbounds F_dist[6 * idx_j - 1] += p_g[11]
        @inbounds F_dist[6 * idx_j    ] += p_g[12]
    end

    # Return the global load vector:
    return F_dist
end