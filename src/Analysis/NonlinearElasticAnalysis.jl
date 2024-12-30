# Reference: 
# S. E. Leon, G. H. Paulino, A. Pereira, I. F. M. Menezes, and E. N. Lages, 
# “A Unified Library of Nonlinear Solution Schemes,” 
# Applied Mechanics Reviews, vol. 64, no. 4, p. 040803, Jul. 2011, DOI: 10.1115/1.4006992.

struct NonlinearElasticAnalysis <: AbstractAnalysisType
    nonlinearsolver::AbstractNonlinearSolver
    maxnumi::Int
    maxnumj::Int
    ϵ::Real
end

struct NonlinearElasticAnalysisCache{
    UT <: Real,
    RT <: Real} <: AbstractSolutionCache
    U::AbstractVector{UT}
    R::AbstractVector{RT}
    states::AbstractVector{ElementState}
end

function solve(model::Model, analysis::NonlinearElasticAnalysis, partitionindices::Vector{Bool})::NonlinearElasticAnalysisCache
    # Initialize the current state of each element:
    elementstates = ElementState[]
    for element in model.elements
        ID = element.ID
        x_i, y_i, z_i = element.node_i.x, element.node_i.y, element.node_i.z
        x_j, y_j, z_j = element.node_j.x, element.node_j.y, element.node_j.z

        L = sqrt((x_j - x_i) ^ 2 + (y_j - y_i) ^ 2 + (z_j - z_i) ^ 2)

        ω = element.ω

        Γ = compute_Γ(x_i, y_i, z_i, x_j, y_j, z_j, ω, ω)

        T     = gettype(element)
        q     = zeros(T, 12)
        k_e_l = zeros(T, 12, 12)
        k_g_l = zeros(T, 12, 12)
        compute_k_e_l!(k_e_l, element, L)
        compute_k_g_l!(k_g_l, element, L, q[7])

        releases_i = element.releases_i
        releases_j = element.releases_j
        condense!(k_e_l, releases_i, releases_j)
        condense!(k_g_l, releases_i, releases_j)

        k_e_g = transform(k_e_l, Γ)
        k_g_g = transform(k_g_l, Γ)
        
        push!(elementstates, ElementState(
            ID,
            x_i, y_i, z_i,
            x_j, y_j, z_j,
            L, 
            ω, ω,
            Γ,
            q,
            k_e_l, k_e_g,
            k_g_l, k_g_g))
    end

    # Assemble the initial global elastic stiffness matrix:
    K_e = assemble_K_e(model, elementstates)
    K_e_ff = K_e[partitionindices, partitionindices]

    # Assemble the global force vector due to concentrated loads:
    F_c   = assemble_F_c(model)
    F_c_f = F_c[partitionindices]

    # Assemble the global force vector due to distributed loads:
    F_d   = assemble_F_d(model)
    F_d_f = F_d[partitionindices]

    # Assemble the reference global force vector:
    P̄ = F_c_f - F_d_f

    # Infer the type of the displacement vector:
    T = promote_type(eltype(K_e), eltype(F_c), eltype(F_d))

    # Preallocate:
    K_e  = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))
    K_g  = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))
    K_t  = zeros(T, count(partitionindices), count(partitionindices))
    δu_p = zeros(T, count(partitionindices))
    δu_r = zeros(T, count(partitionindices))
    δU   = zeros(T, 6 * length(model.nodes))
    U    = zeros(T, 6 * length(model.nodes))
    U_f  = zeros(T, count(partitionindices))
    P_f  = zeros(T, count(partitionindices))
    R_f  = zeros(T, count(partitionindices))
    Q    = zeros(T, 6 * length(model.nodes))
    λ    = zero(T)

    # Extract the nonlinear solver and its parameters:
    nonlinearsolver = analysis.nonlinearsolver
    maxnumi         = analysis.maxnumi
    maxnumj         = analysis.maxnumj
    ϵ               = analysis.ϵ

    # Initialize the increment counter:
    i = 1

    # Enter the increment loop:
    while i ≤ maxnumi
        # Initialize the iteration counter:
        j = 1

        # Initialize the convergance flag:
        converganceflag = false

        # Enter the iteration loop:
        while j ≤ maxnumj && converganceflag ≠ true
            if j == 1
                # Assemble the global elastic stiffness matrix:
                K_e .= 0
                assemble_K_e!(K_e, model, elementstates)

                # Assemble the global geometric stiffness matrix:
                K_g .= 0
                assemble_K_g!(K_g, model, elementstates)

                # Assemble the global tangent stiffness matrix and partition it:
                K_t .= (K_e + K_g)[partitionindices, partitionindices]

                # Compute the displacement increment vector due to P̄ for the free DOFs:
                δu_p .= K_e[partitionindices, partitionindices] \ P̄
            end

            # Compute the displacement increment vector due to R for the free DOFs:
            δu_r .= j == 1 ? 0 : K_t \ R_f
            
            # Compute the load factor increment:
            a = zeros(T, count(partitionindices))
            b = 1
            c = j == 1 ? nonlinearsolver.Δλ : 0
            δλ = constraintequation(a, b, c, δu_p, δu_r)

            # Update the total load factor:
            λ += δλ

            # Update the force vector for the free DOFs:
            P_f += δλ * P̄

            # Update the displacement vector for the free DOFs:
            δU_f = δλ * δu_p + δu_r
            U_f += δU_f

            # Assemble the global displacement vector:
            δU[partitionindices] .= δU_f
            U[partitionindices] .= U_f

            # Update the state of each element:
            for (element, elementstate) in zip(model.elements, elementstates)
                # Extract the nodal displacements increments in the local coordinate system of an element:
                δu_i_g = getnodaldisplacements(model, δU, element.node_i.ID)
                δu_j_g = getnodaldisplacements(model, δU, element.node_j.ID)

                # Update the nodal coordinates of the nodes (i) and (j) of the element:
                elementstate.x_i += δu_i_g[1]
                elementstate.y_i += δu_i_g[2]
                elementstate.z_i += δu_i_g[3]
                elementstate.x_j += δu_j_g[1]
                elementstate.y_j += δu_j_g[2]
                elementstate.z_j += δu_j_g[3]

                elementstate.ω_i += δu_i_g[4]
                elementstate.ω_j += δu_j_g[4]

                Γ′ = compute_Γ(
                    elementstate.x_i, elementstate.y_i, elementstate.z_i, 
                    elementstate.x_j, elementstate.y_j, elementstate.z_j, 
                    elementstate.ω_i,
                    elementstate.ω_j)

                δu_g = [δu_i_g; δu_j_g]
                δu_l = elementstate.Γ * δu_g

                δq_l  = (elementstate.k_e_l + elementstate.k_g_l) * δu_l
                q_l   = elementstate.q + δq_l
                q_g   = transpose(elementstate.Γ) * q_l
                q_l′  = Γ′ * q_g
                elementstate.q = q_l′

                elementstate.Γ = Γ′

                L = sqrt(
                    (elementstate.x_j - elementstate.x_i) ^ 2 + 
                    (elementstate.y_j - elementstate.y_i) ^ 2 + 
                    (elementstate.z_j - elementstate.z_i) ^ 2)
    
                k_e_l = zeros(T, 12, 12)
                k_g_l = zeros(T, 12, 12)
                compute_k_e_l!(k_e_l, element, L)
                compute_k_g_l!(k_g_l, element, L, elementstate.q[7])
    
                releases_i = element.releases_i
                releases_j = element.releases_j
                condense!(k_e_l, releases_i, releases_j)
                condense!(k_g_l, releases_i, releases_j)
    
                elementstate.k_e_l = k_e_l
                elementstate.k_g_l = k_g_l
    
                elementstate.k_e_g = transform(k_e_l, Γ′)
                elementstate.k_g_g = transform(k_g_l, Γ′)
            end

            # Compute the internal force vector:
            Q .= 0
            for (element, elementstate) in zip(model.elements, elementstates)
                index_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
                index_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)

                q = transpose(elementstate.Γ) * elementstate.q
                
                Q[(6 * index_i - 5):(6 * index_i)] += q[ 1:6]
                Q[(6 * index_j - 5):(6 * index_j)] += q[7:12]
            end

            # 
            Q_f = Q[partitionindices]

            # Compute the residual force vector for the free DOFs:
            R_f .= P_f - Q_f

            # Check for convergance:
            converganceflag = norm(R_f) / norm(P̄) ≤ ϵ ? true : false
            
            # Update the iteration counter:
            j += 1
        end

        if converganceflag ≠ true
            @warn "Increment $(i) did not converge."
        end

        # Update the increment counter:
        i += 1
    end

    # Return the solution cache:
    return NonlinearElasticAnalysisCache(U, U, elementstates)
end

