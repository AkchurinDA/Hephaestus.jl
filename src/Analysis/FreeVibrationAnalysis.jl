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
struct FreeVibrationAnalysisCache{
    ΩT <: Real,
    ΦT <: Real} <: AbstractAnalysisCache
    Ω::AbstractVector{ΩT}
    Φ::AbstractMatrix{ΦT}
end

function solve(model::Model, analysistype::FreeVibrationAnalysis, indices_f::Vector{Bool}, indices_s::Vector{Bool})
    # Assemble the global elastic stiffness matrix and partition it:
    K_e    = assemble_K_e(model)
    K_e_ff = K_e[indices_f, indices_f]

    # Assemble the global mass matrix and partition it:
    M    = assemble_M(model)
    M_ff = M[indices_f, indices_f]

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
    return FreeVibrationAnalysisCache(Ω, Φ)
end