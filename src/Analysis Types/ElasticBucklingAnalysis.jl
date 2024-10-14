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
struct ElasticBucklingAnalysisCache <: AbstractAnalysisCache
    Λ            ::AbstractVector{<:Real}
    Φ            ::AbstractMatrix{<:Real}
    K_e          ::AbstractMatrix{<:Real}
    K_g          ::AbstractMatrix{<:Real}
    mapping      ::Dict{Int, Int}
    indices_f    ::Vector{Bool}
    indices_s    ::Vector{Bool}
end

function solve(model::Model, ::ElasticBucklingAnalysis)
    # Perform the linear elastic analysis:
    solution = solve(model, LinearElasticAnalysis())

    # Extract the internal node tags and partitioning indices:
    mapping   = solution.mapping
    indices_f = solution.indices_f
    indices_s = solution.indices_s

    # Extract the global elastic stiffness matrix and partition it:
    K_e = solution.K_e
    K_e_ff = K_e[indices_f, indices_f]
    # K_e_fs = K_e[indices_f, indices_s]
    # K_e_sf = K_e[indices_s, indices_f]
    # K_e_ss = K_e[indices_s, indices_s]
    
    # Compute the internal axial forces in the elements:
    P = [get_element_f_l(model, solution, tag)[7] for (tag, element) in model.elements]

    # Assemble the global geometric stiffness matrix and partition it:
    K_g = _assemble_K_g(model, mapping, P)
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
        Φ[indices_f, i] .= ϕ[:, i] / maximum(abs.(ϕ[:, i]))
    end

    # Return the analysis cache:
    return ElasticBucklingAnalysisCache(Λ, Φ, K_e, K_g, mapping, indices_f, indices_s)
end