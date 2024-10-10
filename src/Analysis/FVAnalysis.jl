"""
    struct FVAnalysis

A type that represents the free vibration (`FV`) analysis.
"""
struct FVAnalysis <: AbstractAnalysisType

end

"""
    struct FVSolutionCache

A type that stores the results of the free vibration analysis.
"""
struct FVSolutionCache <: AbstractSolutionCache
    Ω::AbstractVector{<:Real}
    Φ::AbstractMatrix{<:Real}
end

function solve(model::Model, analysis::FVAnalysis)::FVSolutionCache
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
    # K_e_fs = K_e[indices_f, indices_s]
    # K_e_sf = K_e[indices_s, indices_f]
    # K_e_ss = K_e[indices_s, indices_s]

    # Assemble the global mass matrix and partition it:
    M    = _assemble_M(model, internal_node_IDs)
    M_ff = M[indices_f, indices_f]
    # M_fs = M[indices_f, indices_s]
    # M_sf = M[indices_s, indices_f]
    # M_ss = M[indices_s, indices_s]

    # Solve the generalized eigenvalue problem:
    Ω², ϕ = eigen(K_e_ff, M_ff)

    # Compute the natural frequencies and take only the real parts:
    Ω = sqrt.(real.(Ω²))
    ϕ = real.(ϕ)

    # Normalize the eigenvectors:
    Φ = zeros(eltype(ϕ), 6 * length(model.nodes), length(Ω))
    for i in 1:length(Ω)
        Φ[indices_f, i] = ϕ[:, i] / maximum(abs.(ϕ[:, i]))
    end

    # Return the analysis cache:
    return FVSolutionCache(Ω, Φ)
end