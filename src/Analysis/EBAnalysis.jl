"""
    struct EBAnalysis

A type that represents the elastic buckling (`EB`) analysis.
"""
struct EBAnalysis <: AbstractAnalysisType

end

"""
    struct EBSolutionCache

A type that stores the results of the elastic buckling analysis.
"""
struct EBSolutionCache <: AbstractSolutionCache
    Λ::AbstractVector{<:Real}
    Φ::AbstractMatrix{<:Real}
end

# TODO: 
# Transform the generalized eigenvalue problem into a standard eigenvalue problem. 
# This will allow to use a broader range of eigenvalue solvers and potentially improve the performance.
function solve(model::Model, analysis::EBAnalysis)::EBSolutionCache
    # Perform the 1st-order elastic analysis:
    solution = solve(model, O1EAnalysis())

    # Extract the internal node IDs and partitioning indices:
    internal_node_IDs = solution.internal_node_IDs
    indices_f = solution.indices_f
    indices_s = solution.indices_s

    # Extract the global elastic stiffness matrix and partition it:
    K_e_ff = solution.K_e[indices_f, indices_f]
    # K_e_fs = solution.K_e[indices_f, indices_s]
    # K_e_sf = solution.K_e[indices_s, indices_f]
    # K_e_ss = solution.K_e[indices_s, indices_s]
    
    # Compute the internal axial forces in the elements:
    P = [get_element_f_l(model, solution, element.ID)[7] for element in values(model.elements)]

    # Assemble the global geometric stiffness matrix and partition it:
    K_g = _assemble_K_g(model, solution.internal_node_IDs, P)
    K_g_ff = K_g[indices_f, indices_f]
    # K_g_fs = K_g[indices_f, indices_s]
    # K_g_sf = K_g[indices_s, indices_f]
    # K_g_ss = K_g[indices_s, indices_s]

    # Solve the generalized eigenvalue problem:
    Λ, ϕ = eigen(K_e_ff, -K_g_ff)

    # Take only the real parts:
    Λ = real.(Λ)
    ϕ = real.(ϕ)

    # Normalize the eigenvectors:
    Φ = zeros(eltype(ϕ), 6 * length(model.nodes), length(Λ))
    for i in 1:length(Λ)
        Φ[indices_f, i] = ϕ[:, i] / maximum(abs.(ϕ[:, i]))
    end

    # Return the analysis cache:
    return EBSolutionCache(Λ, Φ)
end