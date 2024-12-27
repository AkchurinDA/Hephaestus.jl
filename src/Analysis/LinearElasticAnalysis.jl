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

function solve(model::Model, analysis::LinearElasticAnalysis, partitionindices::Vector{Bool})
    # Assemble the global geometric stiffness matrix:
    K_e = assemble_K_e(model)
    K_e_ff = K_e[  partitionindices, partitionindices]
    K_e_sf = K_e[.!partitionindices, partitionindices]

    # Assemble the global force vector due to concentrated loads:
    F_conc = assemble_F_conc(model)
    F_conc_f = F_conc[partitionindices]

    # Assemble the global force vector due to distributed loads:
    F_dist = assemble_F_dist(model)
    F_dist_f = F_dist[  partitionindices]
    F_dist_s = F_dist[.!partitionindices]

    # Compute the displacements at the free DOFs:
    U_f = K_e_ff \ (F_conc_f - F_dist_f)

    # Compute the reaction forces at the support DOFs:
    R_s = K_e_sf * U_f + F_dist_s

    # Assemble the global displacement and reaction force vectors:
    U = zeros(eltype(U_f), 6 * length(model.nodes))
    R = zeros(eltype(R_s), 6 * length(model.nodes))
    U[  partitionindices] = U_f
    R[.!partitionindices] = R_s

    # Return the solution cache:
    return LinearElasticAnalysisCache(U, R)
end