struct LinearElasticAnalysis <: AbstractAnalysisType
end

function solve(model::Model, analysis::LinearElasticAnalysis)
    # Extract the partition indices:
    indices_f, indices_s = getpartitionindices(model)

    # Assemble the global geometric stiffness matrix:
    K_e = assemble_K_e(model)
    K_e_ff = K_e[indices_f, indices_f]

    # Assemble the global force vector:
    F = assemble_F(model)
    F_f = F[indices_f]

    # Compute the displacements at the free DOFs:
    U_f = K_e_ff \ F_f

    # Assemble the global displacement vector:
    U = zeros(eltype(U_f), 6 * length(model.nodes))
    U[indices_f] = U_f

    return U
end