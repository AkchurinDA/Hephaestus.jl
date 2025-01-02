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

function solve!(model::Model, analysis::LinearElasticAnalysis, partitionindices::Vector{Bool})::Model
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
    δU_f = K_e_ff \ (F_c_f - F_d_f)

    # Infer type of the displacement vector:
    T = promote_type(eltype(K_e_ff), eltype(F_c), eltype(F_d))

    # TODO: Update the node and element states using the computed displacements.

    # Assemble the global displacement and reaction force vectors:
    δU = zeros(T, 6 * length(model.nodes))
    δU[partitionindices] = δU_f

    # Update the state of the nodes:
    for node in model.nodes
        δu = getnodaldisplacements(model, δU, node.ID)

        updatestate!(node, δu)
    end

    # Update the state of each element:
    for element in model.elements
        δu_i = getnodaldisplacements(model, δU, element.node_i.ID)
        δu_j = getnodaldisplacements(model, δU, element.node_j.ID)

        updatestate!(element, δu_i, δu_j)
    end

    # Return the solution cache:
    return model
end
