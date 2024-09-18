# First-order elastic analysis:
struct O1E <: AbstractAnalysisType

end

struct O1EOutput
    K_e     ::Matrix{<:Real}
    F       ::Vector{<:Real}
    U       ::Vector{<:Real}
    R       ::Vector{<:Real}
end

function solve(model::Model, analysis_type::O1E)
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

    # Assemble the global force vector:
    F = _assemble_F(model, mapping)

    # Partition the global elastic stiffness matrix:
    K_e_ff = K_e[indices_f, indices_f]
    K_e_fs = K_e[indices_f, indices_s]
    K_e_sf = K_e[indices_s, indices_f]
    K_e_ss = K_e[indices_s, indices_s]

    # Partition the global force vector:
    F_f = F[indices_f]

    # Solve for the displacements of the free DOFs:
    if det(K_e_ff) == 0
        throw(ErrorException("The global elastic stiffness matrix is singular. Aborting the analysis."))
    else
        U_f = K_e_ff \ F_f
    end

    # Unpartition the global displacement vector:
    U = _unpartition_U(model, indices_f, U_f)

    # Compute the reaction force vector at the supported DOFs:
    R = K_e_sf * U_f

    # Return the displacements of the free DOFs and the reaction force vector:
    return O1EOutput(K_e, F, R, U)
end
