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
    ΦT <: Real} <: AbstractSolutionCache
    Λ::AbstractVector{ΛT}
    Φ::AbstractMatrix{ΦT}
end

function solve(model::Model, ::ElasticBucklingAnalysis, partitionindices::Vector{Bool})
    # Perform the linear elastic analysis:
    solution = solve(model, LinearElasticAnalysis(), partitionindices)

    # Extract the global elastic stiffness matrix and partition it:
    K_e = assemble_K_e(model)
    K_e_ff = K_e[partitionindices, partitionindices]

    # Extract the element axial loads and partition them:
    P = [getelementforces(model, solution, element.ID)[7] for element in model.elements]

    # Assemble the global geometric stiffness matrix and partition it:
    K_g = assemble_K_g(model, P)
    K_g_ff = K_g[partitionindices, partitionindices]

    # Solve the generalized eigenvalue problem:
    Λ, ϕ = eigen(K_e_ff, -K_g_ff)

    # Take only the real parts:
    Λ = real.(Λ)
    ϕ = real.(ϕ)

    # Normalize the eigenvectors:
    Φ = zeros(eltype(ϕ), 6 * length(model.nodes), length(Λ))
    for i in 1:length(Λ)
        Φ[partitionindices, i] .= ϕ[:, i] / maximum(abs.(ϕ[:, i]))
    end

    # Return the analysis cache:
    return ElasticBucklingAnalysisCache(Λ, Φ)
end