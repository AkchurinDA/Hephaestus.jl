struct EB <: AbstractAnalysisType
    
end

struct EBCache
    K_e
    K_g
    λ
    Φ
end

function solve(model::Model, analysis_type::EB)
    # Assign internal tags to nodes:
    mapping = Dict{Int, Int}()
    for (i, node) in enumerate(values(model.nodes))
        mapping[node.tag] = i
    end

    # Assign internal tags to elements:
    element_tag_mapping = Dict{Int, Int}()
    for (i, element) in enumerate(values(model.elements))
        element_tag_mapping[element.tag] = i
    end

    # Get the partition indices:
    indices_f, indices_s = _get_partition_indices(model, mapping)

    # Assemble the global elastic stiffness matrix:
    K_e = _assemble_K_e(model, mapping)

    # Assemble the global geometric stiffness matrix:
    K_g = _assemble_K_g(model, mapping)

    # Partition the global elastic stiffness matrix:
    K_e_ff = K_e[indices_f, indices_f]

    # Partition the global geometric stiffness matrix:
    K_g_ff = K_g[indices_f, indices_f]

    # Solve the eigenvalue problem:
    λ, Φ_f = eigen(K_e_ff, K_g_ff)

    # Normalize the eigenvectors:
    for i in axes(Φ_f, 2)
        Φ_f[:, i] /= maximum(abs.(Φ_f[:, i]))
    end

    # Unpartition the global displacement vector using the eigenvectors:
    Φ = zeros(eltype(Φ_f), 6 * length(model.nodes), size(Φ_f, 2))
    for i in axes(Φ_f, 2)
        Φ[:, i] = _unpartition_U(model, indices_f, Φ_f[:, i])
    end

    return EBCache(K_e, K_g, λ, Φ)
end