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

function solve(model::Model, analysis::LinearElasticAnalysis, indices_f::Vector{Bool}, indices_s::Vector{Bool})
    # Extract the partition indices:
    indices_f, indices_s = getpartitionindices(model)

    # Assemble the global geometric stiffness matrix:
    K_e = assemble_K_e(model)
    K_e_ff = K_e[indices_f, indices_f]
    K_e_sf = K_e[indices_s, indices_f]

    # Assemble the global force vector due to concentrated loads:
    F_conc = assemble_F_conc(model)
    F_conc_f = F_conc[indices_f]

    # Assemble the global force vector due to distributed loads:
    F_dist = assemble_F_dist(model)
    F_dist_f = F_dist[indices_f]
    F_dist_s = F_dist[indices_s]

    # Compute the displacements at the free DOFs:
    U_f = K_e_ff \ (F_conc_f - F_dist_f)

    # Assemble the global displacement vector:
    U = zeros(eltype(U_f), 6 * length(model.nodes))
    U[indices_f] = U_f

    # Compute the reaction forces at the support DOFs:
    R_s = K_e_sf * U_f + F_dist_s

    # Assemble the global reaction force vector:
    R = zeros(eltype(R_s), 6 * length(model.nodes))
    R[indices_s] = R_s

    # Return the solution cache:
    return LinearElasticAnalysisCache(U, R)
end