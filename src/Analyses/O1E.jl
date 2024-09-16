# First-order elastic analysis:
struct O1E <: AbstractAnalysisType

end

struct O1EOutput
    K_e     ::Matrix{<:Real}
    K_e_ff  ::Matrix{<:Real}
    K_e_fs  ::Matrix{<:Real}
    K_e_sf  ::Matrix{<:Real}
    K_e_ss  ::Matrix{<:Real}
    F       ::Vector{<:Real}
    F_f     ::Vector{<:Real}
    F_s     ::Vector{<:Real}
    U       ::Vector{<:Real}
    U_f     ::Vector{<:Real}
    U_s     ::Vector{<:Real}
end

function solve(model::Model, analysis_type::O1E)
    # Assign internal tags to nodes:
    node_tag_mapping = Dict{Int, Int}()
    for (i, node) in enumerate(values(model.nodes))
        node_tag_mapping[node.tag] = i
    end

    # Assign internal tags to elements:
    element_tag_mapping = Dict{Int, Int}()
    for (i, element) in enumerate(values(model.elements))
        element_tag_mapping[element.tag] = i
    end

    # Get the partition indices:
    indices_f, indices_s = _get_partition_indices(model, node_tag_mapping)
    
    # Assemble the global elastic stiffness matrix:
    K_e = _assemble_K_e(model, node_tag_mapping)

    # Assemble the global force vector:
    F = _assemble_F(model, node_tag_mapping)

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

    # Compute the reaction force vector at the supported DOFs:
    F_s = K_e_sf * U_f

    # Return the displacements of the free DOFs and the reaction force vector:
    return O1EOutput(K_e, K_e_ff, K_e_fs, K_e_sf, K_e_ss, F, F_f, F_s, U_f, U_f, U_f)
end

function _get_partition_indices(model::Model, node_tag_mapping::Dict{Int, Int})
    # Preallocate:
    indices_f = Int[] # List of free DOFs
    indices_s = Int[] # List of supported DOFs

    # Get the partition indices:
    for node in values(model.nodes)
        node_tag_i = node_tag_mapping[node.tag]

        if haskey(model.supports, node.tag) # Check if the node's DOFs are fixed
            for i in 1:6
                model.supports[node.tag][i] ? push!(indices_s, 6 * (node_tag_i - 1) + i) : push!(indices_f, 6 * (node_tag_i - 1) + i)
            end
        else
            for i in 1:6
                push!(indices_f, 6 * (node_tag_i - 1) + i)
            end
        end
    end

    # Return the partition indices:
    return indices_f, indices_s
end

function _assemble_K_e(model::Model, node_tag_mapping::Dict{Int, Int})
    # Preallocate the global elastic stiffness matrix:
    T   = promote_type([eltype(element.k_e_g) for element in values(model.elements)]...)
    K_e = zeros(T, 6 * length(model.nodes), 6 * length(model.nodes))

    # Assemble the global elastic stiffness matrix:
    for element in values(model.elements)
        # Extract the internal node IDs of the nodes of an element:
        node_i_tag_i = node_tag_mapping[element.node_i_tag]
        node_j_tag_i = node_tag_mapping[element.node_j_tag]
 
        range_i = (6 * (node_i_tag_i - 1) + 1):(6 * node_i_tag_i)
        range_j = (6 * (node_j_tag_i - 1) + 1):(6 * node_j_tag_i)

        @inbounds K_e[range_i, range_i] += element.k_e_g[1:6 , 1:6 ]
        @inbounds K_e[range_i, range_j] += element.k_e_g[1:6 , 7:12]
        @inbounds K_e[range_j, range_i] += element.k_e_g[7:12, 1:6 ]
        @inbounds K_e[range_j, range_j] += element.k_e_g[7:12, 7:12]
    end

    # Return the global elastic stiffness matrix:
    return K_e
end

function _assemble_F(model::Model, node_tag_mapping::Dict{Int, Int})
    # Preallocate the global force vector:
    T = promote_type([eltype(nodal_loads) for nodal_loads in values(model.nodal_loads)]...)
    F = zeros(T, 6 * length(model.nodes))

    # Assemble the global force vector:
    for node in values(model.nodes)
        if haskey(model.nodal_loads, node.tag)
            node_tag_i = node_tag_mapping[node.tag]

            for i in 1:6
                if model.nodal_loads[node.tag][i] != 0
                    @inbounds F[6 * (node_tag_i - 1) + i] = model.nodal_loads[node.tag][i]
                end
            end
        end
    end

    # Return the global force vector:
    return F
end