# Reference:
# S. E. Leon, G. H. Paulino, A. Pereira, I. F. M. Menezes, and E. N. Lages,
# “A Unified Library of Nonlinear Solution Schemes,”
# Applied Mechanics Reviews, vol. 64, no. 4, p. 040803, Jul. 2011, DOI: 10.1115/1.4006992.

"""
    struct NonlinearElasticAnalysis <: AbstractAnalysisType

A type representing (geometrically) nonlinear (materially) elastic analysis.
"""
struct NonlinearElasticAnalysis <: AbstractAnalysisType
    nonlinearsolver::AbstractNonlinearSolver
    maxnumi::Int
    maxnumj::Int
    update::Symbol
    ϵ::Real
end

function solve!(model::Model, analysis::NonlinearElasticAnalysis, partitionindices::Vector{Bool})::Model
    # Infer the type of the global elastic stiffness matrix:
    KT = promote_type([gettype(element) for element in model.elements]...)

    # Assemble the global force vector due to concentrated loads:
    F_c   = assemble_F_c(model)
    F_c_f = F_c[partitionindices]

    # Assemble the global force vector due to distributed loads:
    F_d   = assemble_F_d(model)
    F_d_f = F_d[  partitionindices]
    F_d_s = F_d[.!partitionindices]

    # Assemble the reference global force vector:
    P̄ = F_c_f - F_d_f

    # Infer the type of the displacement vector:
    T = promote_type(KT, eltype(F_c), eltype(F_d))

    # Preallocate:
    K_e  = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))
    K_g  = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))
    K_t_ff = zeros(T, count(  partitionindices), count(partitionindices))
    K_t_sf = zeros(T, count(.!partitionindices), count(partitionindices))
    δu_p = zeros(T, count(partitionindices))
    δu_r = zeros(T, count(partitionindices))
    δU   = zeros(T, 6 * length(model.nodes))
    δR   = zeros(T, 6 * length(model.nodes))
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
                # Assemble the global elastic stiffness matrix:
                K_e .= 0
                assemble_K_e!(K_e, model)

                # Assemble the global geometric stiffness matrix:
                K_g .= 0
                assemble_K_g!(K_g, model)

                # Assemble the global tangent stiffness matrix and partition it:
                K_t_ff .= K_e[  partitionindices, partitionindices] + K_g[  partitionindices, partitionindices]
                K_t_sf .= K_e[.!partitionindices, partitionindices] + K_g[.!partitionindices, partitionindices]

                # Compute the displacement increment vector due to P̄ for the free DOFs:
                δu_p .= K_t_ff \ P̄
            end

            # Compute the displacement increment vector due to R for the free DOFs:
            δu_r .= j == 1 ? 0 : K_t_ff \ R_f

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

            # Comptute the global reaction force vector:
            δR_s = K_t_sf * δU_f + δλ * F_d_s

            # Assemble the global reaction vector:
            δR[.!partitionindices] .= δR_s

            # Update the state of the nodes:
            for node in model.nodes
                δu = getnodaldisplacements(model, δU, node.ID)
                δr = getnodalreactions(model, δR, node.ID)

                updatestate!(node, δu, δr)
            end

            # Update the state of each element:
            for element in model.elements
                δu_i = getnodaldisplacements(model, δU, element.node_i.ID)
                δu_j = getnodaldisplacements(model, δU, element.node_j.ID)

                updatestate!(element, δu_i, δu_j)
            end

            # Compute the global internal force vector:
            Q .= 0
            for element in model.elements
                index_i = findfirst(x -> x.ID == element.node_i.ID, model.nodes)
                index_j = findfirst(x -> x.ID == element.node_j.ID, model.nodes)

                q = transpose(element.state.Γ) * element.state.q

                Q[(6 * index_i - 5):(6 * index_i)] += q[ 1:6]
                Q[(6 * index_j - 5):(6 * index_j)] += q[7:12]
            end

            # Partition the global internal force vector:
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
    return model
end
