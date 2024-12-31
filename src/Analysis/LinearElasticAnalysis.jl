"""
    struct LinearElasticAnalysis

A type representing the (geometrically) linear (materially) elastic analysis.
"""
struct LinearElasticAnalysis <: AbstractAnalysisType

end

"""
    struct LinearElasticAnalysisCache

A type used to store the results of linear elastic analysis.
"""
struct LinearElasticAnalysisCache{
    UT <: Real,
    RT <: Real} <: AbstractSolutionCache
    U::AbstractVector{UT}
    R::AbstractVector{RT}
end

function solve(model::Model, analysis::LinearElasticAnalysis, partitionindices::Vector{Bool})::LinearElasticAnalysisCache
    # Assemble the global geometric stiffness matrix:
    K_e    = assemble_K_e(model)
    K_e_ff = K_e[  partitionindices, partitionindices]
    K_e_sf = K_e[.!partitionindices, partitionindices]

    # Check if the global stiffness matrix is singular:
    if det(K_e_ff) ≈ 0
        error("The global stiffness matrix is singular. Aborting.")
    end
    
    # Assemble the global force vector due to concentrated loads:
    F_c = assemble_F_c(model)
    F_c_f = F_c[partitionindices]

    # Assemble the global force vector due to distributed loads:
    F_d   = assemble_F_d(model)
    F_d_f = F_d[  partitionindices]
    F_d_s = F_d[.!partitionindices]

    # Compute the displacements at the free DOFs:
    U_f = K_e_ff \ (F_c_f - F_d_f)    

    # TODO: Update the node and element states using the computed displacements.

    # Compute the reaction forces at the support DOFs:
    R_s = K_e_sf * U_f + F_d_s

    # Assemble the global displacement and reaction force vectors:
    U = zeros(eltype(U_f), 6 * length(model.nodes))
    R = zeros(eltype(R_s), 6 * length(model.nodes))
    U[  partitionindices] = U_f
    R[.!partitionindices] = R_s

    # Return the solution cache:
    return LinearElasticAnalysisCache(U, R)
end