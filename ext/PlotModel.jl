Makie.@recipe(PlotModel, Model) do scene
    Makie.Attributes(
        # Nodes:
        show_nodes       = true,
        node_color       = :steelblue,
        node_marker      = :circle,
        node_markersize  = 6,
        node_strokecolor = :black,
        node_strokewidth = 1,

        # Node labels:
        show_node_labels = true,
        node_label_align = (:center, :bottom),
        node_label_color = :steelblue,

        # Elements:
        show_elements     = true,
        element_color     = :black,
        element_linestyle = :solid,
        element_linewidth = 1,
        
        # Element labels:
        show_element_labels = true,
        element_label_align = (:center, :bottom),
        element_label_color = :black,
        
        # Supports:
        show_supports = false,
        support_color = :green)
end

Makie.preferred_axis_type(::PlotModel) = Makie.Axis3

function Makie.plot!(P::PlotModel)
    model = P[:Model]

    if !isempty(model[].elements)
        if P[:show_elements][]
            for element in values(model[].elements)
                x_i, y_i, z_i = element.x_i, element.y_i, element.z_i
                x_j, y_j, z_j = element.x_j, element.y_j, element.z_j 

                lines!(P, [x_i, x_j], [y_i, y_j], [z_i, z_j], 
                    color     = P[:element_color    ],
                    linestyle = P[:element_linestyle],
                    linewidth = P[:element_linewidth])

                if P[:show_element_labels][]
                    x_m = (x_i + x_j) / 2
                    y_m = (y_i + y_j) / 2
                    z_m = (z_i + z_j) / 2

                    text!(P, [x_m], [y_m], [z_m],
                        text  = string(element.ID),
                        color = P[:element_label_color],
                        align = P[:element_label_align])
                end
            end
        end
    else
        @warn "The model has no elements to plot."
    end

    if !isempty(model[].nodes)
        if P[:show_nodes][]
            for node in values(model[].nodes)
                x, y, z = node.x, node.y, node.z

                scatter!(P, [x], [y], [z], 
                    color       = P[:node_color      ],
                    marker      = P[:node_marker     ],
                    markersize  = P[:node_markersize ],
                    strokecolor = P[:node_strokecolor],
                    strokewidth = P[:node_strokewidth],
                    overdraw    = true)

                if P[:show_node_labels][]
                    text!(P, [x], [y], [z],
                        text  = string(node.ID),
                        color = P[:node_label_color],
                        align = P[:node_label_align])
                end
            end
        end
    else
        @warn "The model has no nodes to plot."
    end

    # Return the updated model plot:
    return P
end