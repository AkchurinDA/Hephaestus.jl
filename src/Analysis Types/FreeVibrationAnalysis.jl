"""
    struct LinearElasticAnalysis

A type representing the free vibration analysis.
"""
struct FreeVibrationAnalysis <: AbstractAnalysisType

end

"""
    struct LinearElasticAnalysisCache

A type used to store the results of free vibration analysis.
"""
struct FreeVibrationAnalysisCache <: AbstractAnalysisCache
    Ω            ::AbstractVector{<:Real}
    Φ            ::AbstractMatrix{<:Real}
    K_e          ::AbstractMatrix{<:Real}
    M            ::AbstractMatrix{<:Real}
    mapping      ::Dict{Int, Int}
    indices_f    ::Vector{Bool}
    indices_s    ::Vector{Bool}
end

function solve(model::Model, ::FreeVibrationAnalysis)
    # Assign the internal node tags:
    mapping = Dict{Int, Int}()
    for (i, tag) in enumerate(keys(model.nodes))
        mapping[tag] = i
    end

    # Get the partition indices:
    indices_f, indices_s = _get_partition_indices(model, mapping)

    # Assemble the global elastic stiffness matrix and partition it:
    K_e    = _assemble_K_e(model, mapping)
    K_e_ff = K_e[indices_f, indices_f]
    # K_e_fs = K_e[indices_f, indices_s]
    # K_e_sf = K_e[indices_s, indices_f]
    # K_e_ss = K_e[indices_s, indices_s]

    # Assemble the global mass matrix and partition it:
    M    = _assemble_M(model, mapping)
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
    return FreeVibrationAnalysisCache(Ω, Φ, K_e, M, mapping, indices_f, indices_s)
end