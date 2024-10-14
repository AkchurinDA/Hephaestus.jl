"""
    struct LinearElasticAnalysis

A type representing the (geometrically) linear (materially) elastic analysis.
"""
struct LinearElasticAnalysis <: AbstractAnalysisType

end

"""
    struct LinearElasticAnalysisCache

A type used to store the results of linear elastic analysis.
"""
struct LinearElasticAnalysisCache <: AbstractAnalysisCache
    K_e          ::AbstractMatrix{<:Real}
    F            ::AbstractVector{<:Real}
    P            ::AbstractVector{<:Real}
    U            ::AbstractVector{<:Real}
    R            ::AbstractVector{<:Real}
    mapping      ::Dict{Int, Int}
    indices_f    ::Vector{Bool}
    indices_s    ::Vector{Bool}
end

function solve(model::Model, ::LinearElasticAnalysis)
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
    K_e_sf = K_e[indices_s, indices_f]
    # K_e_ss = K_e[indices_s, indices_s]

    # Assemble the global force vector and partition it:
    F   = _assemble_F(model, mapping)
    F_f = F[indices_f]
    # F_s = F[indices_s]

    # Assemble the global fixed-end force vector and partition it:
    P   = _assemble_P(model, mapping)
    P_f = P[indices_f]
    P_s = P[indices_s]

    # Compute the displacements at the free DOFs:
    U_f = K_e_ff \ (F_f - P_f)

    # Assemble the global displacement vector:
    U             = zeros(eltype(U_f), 6 * length(model.nodes))
    U[indices_f] .= U_f

    # Compute the reaction forces at the supported DOFs:
    R_s = K_e_sf * U_f + P_s

    # Assemble the global reaction force vector:
    R             = zeros(eltype(R_s), 6 * length(model.nodes))
    R[indices_s] .= R_s

    # Return the analysis cache:
    return LinearElasticAnalysisCache(K_e, F, P, U, R, mapping, indices_f, indices_s)
end