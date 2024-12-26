"""
    struct ElasticBucklingAnalysis

A type representing the elastic buckling analysis.
"""
struct ElasticBucklingAnalysis <: AbstractAnalysisType

end

"""
    struct ElasticBucklingAnalysisCache

A type used to store the results of elastic buckling analysis.
"""
struct ElasticBucklingAnalysisCache{
    ΛT <: Real,
    ΦT <: Real} <: AbstractAnalysisCache
    Λ::AbstractVector{ΛT}
    Φ::AbstractMatrix{ΦT}
end

function solve(model::Model, ::ElasticBucklingAnalysis, indices_f::Vector{Bool}, indices_s::Vector{Bool})
    # Perform the linear elastic analysis:
    solution = solve(model, LinearElasticAnalysis(), indices_f, indices_s)

    # Extract the axial loads within each element:
    # TODO

    # Extract the global elastic stiffness matrix and partition it:
    K_e = assemble_K_e(model)
    K_e_ff = K_e[indices_f, indices_f]

    # Assemble the global geometric stiffness matrix and partition it:
    K_g = assemble_K_g(model, P)
    K_g_ff = K_g[indices_f, indices_f]

    # Solve the generalized eigenvalue problem:
    Λ, ϕ = eigen(K_e_ff, -K_g_ff)

    # Take only the real parts:
    Λ = real.(Λ)
    ϕ = real.(ϕ)

    # Normalize the eigenvectors:
    Φ = zeros(eltype(ϕ), 6 * length(model.nodes), length(Λ))
    for i in 1:length(Λ)
        Φ[indices_f, i] .= ϕ[:, i] / maximum(abs.(ϕ[:, i]))
    end

    # Return the analysis cache:
    return ElasticBucklingAnalysisCache(Λ, Φ)
end