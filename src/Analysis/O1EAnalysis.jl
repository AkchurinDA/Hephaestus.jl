"""
    struct O1EAnalysis

A type that represents the 1st-order (`O1`) elastic (`E`) analysis.
"""
struct O1EAnalysis <: AbstractAnalysisType

end

"""
    struct O1ESolutionCache

A type that stores the results of the 1st-order elastic analysis.
"""
struct O1ESolutionCache <: AbstractSolutionCache
    K_e                 ::AbstractMatrix{<:Real}
    F                   ::AbstractVector{<:Real}
    P                   ::AbstractVector{<:Real}
    U                   ::AbstractVector{<:Real}
    R                   ::AbstractVector{<:Real}

    internal_node_IDs   ::Dict{Int, Int}
    indices_f           ::Vector{Int}
    indices_s           ::Vector{Int}
end

"""
    solve(model::Model, analysis::O1EAnalysis)

Performs the 1st-order elastic analysis on the model.
"""
function solve(model::Model, analysis::O1EAnalysis)::O1ESolutionCache
    # Assign internal IDs to the nodes:
    internal_node_IDs = Dict{Int, Int}()
    for (i, node) in enumerate(values(model.nodes))
        internal_node_IDs[node.ID] = i
    end

    # Compute the partitioning indices:
    indices_f, indices_s = _get_partitioning_indices(model, internal_node_IDs)

    # Assemble the global elastic stiffness matrix and partition it:
    K_e    = _assemble_K_e(model, internal_node_IDs)
    K_e_ff = K_e[indices_f, indices_f]
    K_e_fs = K_e[indices_f, indices_s]
    K_e_sf = K_e[indices_s, indices_f]
    K_e_ss = K_e[indices_s, indices_s]

    # Assemble the global force vector and partition it:
    F = _assemble_F(model, internal_node_IDs)
    F_f = F[indices_f]
    F_s = F[indices_s]

    # Assemble the global fixed-end force vector and partition it:
    P = _assemble_P(model, internal_node_IDs)
    P_f = P[indices_f]
    P_s = P[indices_s]

    # Compute the displacements at the free DOFs:
    U_f = K_e_ff \ (F_f - P_f)

    # Assemble the global displacement vector:
    U = zeros(eltype(U_f), 6 * length(model.nodes))
    U[indices_f] = U_f

    # Compute the reaction forces at the supported DOFs:
    R_s = K_e_sf * U_f + P_s

    # Assemble the global reaction force vector:
    R = zeros(eltype(R_s), 6 * length(model.nodes))
    R[indices_s] = R_s

    # Return the analysis cache:
    return O1ESolutionCache(K_e, F, P, U, R, internal_node_IDs, indices_f, indices_s)
end