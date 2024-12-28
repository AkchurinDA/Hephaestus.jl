# Reference: 
# S. E. Leon, G. H. Paulino, A. Pereira, I. F. M. Menezes, and E. N. Lages, 
# “A Unified Library of Nonlinear Solution Schemes,” 
# Applied Mechanics Reviews, vol. 64, no. 4, p. 040803, Jul. 2011, DOI: 10.1115/1.4006992.

struct NonlinearElasticAnalysis <: AbstractAnalysisType
    nonlinearsolver::AbstractNonlinearSolver
    maxnumi::Int
    maxnumj::Int
    update::Symbol
    ϵ::Real
end

struct NonlinearElasticAnalysisCache{
    UT <: Real,
    RT <: Real} <: AbstractSolutionCache
    U::AbstractVector{UT}
    R::AbstractVector{RT}
end

function solve(model::Model, analysis::NonlinearElasticAnalysis, partitionindices::Vector{Bool})
    # Assemble the initial global elastic stiffness matrix:
    K_e    = assemble_K_e(model)
    K_e_ff = K_e[partitionindices, partitionindices]

    # Assemble the global force vector due to concentrated loads:
    F_conc   = assemble_F_conc(model)
    F_conc_f = F_conc[partitionindices]

    # Assemble the global force vector due to distributed loads:
    F_dist   = assemble_F_dist(model)
    F_dist_f = F_dist[  partitionindices]

    # Assemble the reference global force vector:
    P̄ = F_conc_f - F_dist_f

    # Infer the type of the displacement vector:
    T = promote_type(eltype(K_e_ff), eltype(P̄))

    # Preallocate:
    K_e  = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))
    K_g  = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))
    K_t  = zeros(T, count(partitionindices), count(partitionindices))
    δu_p = zeros(T, count(partitionindices))
    δu_r = zeros(T, count(partitionindices))
    dU   = zeros(T, 6 * length(model.nodes))
    U    = zeros(T, 6 * length(model.nodes))
    U_f  = zeros(T, count(partitionindices))
    P_f  = zeros(T, count(partitionindices))
    R_f  = zeros(T, count(partitionindices))
    Q    = zeros(T, 6 * length(model.nodes))
    λ    = zero(T)

    # Initialize the current state of each element:
    elementstates = [ElementState(
        element.ID,
        element.Γ, 
        element.ω,
        element.ω,
        [element.node_i.x, element.node_i.y, element.node_i.z, 0, 0, 0],
        [element.node_j.x, element.node_j.y, element.node_j.z, 0, 0, 0],
        0) for element in model.elements]

    # Extract the nonlinear solver and its parameters:
    nonlinearsolver = analysis.nonlinearsolver
    maxnumi         = analysis.maxnumi
    maxnumj         = analysis.maxnumj
    update          = analysis.update
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
            if j == 1 || update == :standard
                # Extract the element axial loads:
                N = [elementstate.N for elementstate in elementstates]
                
                # Assemble the global elastic stiffness matrix:
                K_e .= 0
                assemble_K_e!(K_e, model)

                # Assemble the global geometric stiffness matrix:
                K_g .= 0
                assemble_K_g!(K_g, model, N)

                # Assemble the global tangent stiffness matrix and partition it:
                K_t .= (K_e + K_g)[partitionindices, partitionindices]

                # Compute the displacement increment vector due to P̄ for the free DOFs:
                δu_p .= K_t \ P̄
            end

            # Compute the displacement increment vector due to R for the free DOFs:
            if j == 1
                δu_r .= 0
            else
                δu_r .= K_t \ R_f
            end
            
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
            dU_f = δλ * δu_p + δu_r
            U_f += dU_f

            # Assemble the global displacement vector:
            dU[partitionindices] .= dU_f
            U[partitionindices] .= U_f

            # Update the state of each element:
            for (elementstate, element) in zip(elementstates, model.elements)
                d_node_i_disp = getnodaldisplacements(model, dU, element.node_i.ID)
                d_node_j_disp = getnodaldisplacements(model, dU, element.node_j.ID)
                
                # Compute the transformation matrix from the previous state to the current state:
                Γ′ = compute_Γ(
                    elementstate.node_i_coords[1:3]..., 
                    elementstate.node_j_coords[1:3]...,
                    d_node_i_disp...,
                    d_node_j_disp...,
                    elementstate.ω_i,
                    elementstate.ω_j)

                elementstate.Γ = Γ′ * element.Γ

                # Update the current element orientation:
                elementstate.ω_i += d_node_i_disp[4]
                elementstate.ω_j += d_node_j_disp[4]

                elementstate.node_i_coords += d_node_i_disp
                elementstate.node_j_coords += d_node_j_disp

                elementstate.N = getelementforces(element, elementstate)[7]

                @show d_node_i_disp
                @show d_node_j_disp
                @show elementstate.N
            end

            # Compute the internal force vector:
            # TODO: To implement this, I need to be able to remember the current state each element.

            # Compute the residual force vector for the free DOFs:
            R_f .= P_f - Q[partitionindices]

            # Check for convergance:
            if norm(R_f) / norm(P̄) ≤ ϵ
                converganceflag = true
            end
            
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
    return NonlinearElasticAnalysisCache(U, U)
end

