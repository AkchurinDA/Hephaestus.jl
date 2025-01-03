"""
    struct FreeVibrationAnalysis <: AbstractAnalysisType

A type representing the free vibration analysis.

To perform an free vibration analysis, use the following command:
```julia
solution = solve(model, FreeVibrationAnalysis())
```
"""
@kwdef struct FreeVibrationAnalysis <: AbstractAnalysisType
    nummodes::Int = 10
end

"""
    struct FreeVibrationAnalysisCache <: AbstractSolutionCache

A type used to store the results of free vibration analysis.
"""
struct FreeVibrationAnalysisCache{
    ΩT <: Real,
    ΦT <: Real} <: AbstractSolutionCache
    Ω::AbstractVector{ΩT}
    Φ::AbstractMatrix{ΦT}
end

function solve(model::Model, ::FreeVibrationAnalysis, partitionindices::Vector{Bool})::FreeVibrationAnalysisCache
    # Assemble the global elastic stiffness matrix and partition it:
    K_e    = assemble_K_e(model)
    K_e_ff = K_e[partitionindices, partitionindices]

    # Assemble the global mass matrix and partition it:
    M    = assemble_M(model)
    M_ff = M[partitionindices, partitionindices]

    # Solve the generalized eigenvalue problem:
    Ω², ϕ = eigen(K_e_ff, M_ff)

    # Compute the natural frequencies and take only the real parts:
    Ω = sqrt.(real.(Ω²))
    ϕ = real.(ϕ)

    # Normalize the eigenvectors:
    Φ = zeros(eltype(ϕ), 6 * length(model.nodes), length(Ω))
    for i in 1:length(Ω)
        Φ[partitionindices, i] = ϕ[:, i] / maximum(abs.(ϕ[:, i]))
    end

    # Return the analysis cache:
    return FreeVibrationAnalysisCache(Ω, Φ)
end
