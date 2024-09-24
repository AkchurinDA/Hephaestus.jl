struct O1EAnalysis <: AbstractAnalysisType
    
end

struct O1ECache <: AbstractAnalysisCache
    internal_node_IDs   ::Dict{Int, Int}
    indices_f           ::Vector{Int}
    indices_s           ::Vector{Int}
    K_e                 ::Matrix{<:Real}
    F                   ::Vector{<:Real}
    P                   ::Vector{<:Real}
    U                   ::Vector{<:Real}
    R                   ::Vector{<:Real}
end

function solve(model::Model, analysis::O1EAnalysis, log::Bool = false)
    # Assign internal node IDs to the nodes and their DOFs:
    internal_node_IDs = Dict{Int, Int}()
    internal_DOF_IDs  = Dict{Int, Int}()
    for (i, node_ID) in enumerate(keys(model.nodes))
        internal_node_IDs[node_ID] = i

        user_base_index     = 6 * (node_ID - 1)
        internal_base_index = 6 * (i - 1)

        for i in 1:6
            internal_DOF_IDs[internal_base_index + i] = user_base_index + i
        end
    end

    # Compute the partitioning indices:
    indices_f, indices_s = _get_partition_indices(model, internal_node_IDs)

    # Assemble the global elastic stiffness matrix and partition it:
    K_e    = _assemble_K_e(model, internal_node_IDs)
    K_e_ff = K_e[indices_f, indices_f]
    K_e_fs = K_e[indices_f, indices_s]
    K_e_sf = K_e[indices_s, indices_f]
    K_e_ss = K_e[indices_s, indices_s]

    # Assemble the global force vector:
    F   = _assemble_F(model, internal_node_IDs)
    F_f = F[indices_f]
    F_s = F[indices_s]

    # Assemble the global fixed-end forces vector:
    P   = _assemble_P(model, internal_node_IDs)
    P_f = P[indices_f]
    P_s = P[indices_s]

    # Compute the global displacement vector:
    U_f = K_e_ff \ (F_f - P_f)

    # Assemble the global displacement vector:
    U = zeros(eltype(U_f), 6 * length(model.nodes))
    for (i, index) in enumerate(indices_f)
        U[internal_DOF_IDs[index]] = U_f[i]
    end

    # Compute the global reaction forces vector:
    R_s = K_e_sf * U_f + P_s

    # Assemble the global reaction force vector:
    R = zeros(eltype(R_s), 6 * length(model.nodes))
    for (i, index) in enumerate(indices_s)
        R[internal_DOF_IDs[index]] = R_s[i]
    end

    # Return the solution:
    return O1ECache(internal_node_IDs, indices_f, indices_s, K_e, F, P, U, R)
end

