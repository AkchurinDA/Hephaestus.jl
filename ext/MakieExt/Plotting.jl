using MakieCore

MakieCore.@recipe(PlotModel, Model) do scene
    MakieCore.Attributes(
        # Nodes:
        node_color       = :red ,
        node_marker      = :circle,
        node_markersize  = 6,
        node_strokecolor = :black,
        node_strokewidth = 1,

        # Node labels:
        label_nodes       = true,
        align_node_labels = (:center, :bottom),
        node_label_color  = :red,

        # Elements:
        element_color     = :black,
        element_linestyle = :solid,
        element_linewidth = 1,
        
        # Element labels:
        label_elements       = true,
        align_element_labels = (:center, :bottom),
        element_label_color  = :black)
end

function MakieCore.plot!(P::PlotModel)
    model = P[:Model]

    for element in values(model[].elements)
        x_i, y_i, z_i = element.x_i, element.y_i, element.z_i
        x_j, y_j, z_j = element.x_j, element.y_j, element.z_j 

        lines!(P, [x_i, x_j], [y_i, y_j], [z_i, z_j], 
            color     = P[:element_color    ],
            linestyle = P[:element_linestyle],
            linewidth = P[:element_linewidth])

        if P[:label_elements][]
            x_m = (x_i + x_j) / 2
            y_m = (y_i + y_j) / 2
            z_m = (z_i + z_j) / 2

            text!(P, [x_m], [y_m], [z_m],
                text  = string(element.ID),
                color = P[:element_label_color ],
                align = P[:align_element_labels])
        end
    end

    for node in values(model[].nodes)
        x, y, z = node.x, node.y, node.z

        scatter!(P, [x], [y], [z], 
            color       = P[:node_color      ],
            marker      = P[:node_marker     ],
            markersize  = P[:node_markersize ],
            strokecolor = P[:node_strokecolor],
            strokewidth = P[:node_strokewidth],
            overdraw    = true)

        if P[:label_nodes][]
            text!(P, [x], [y], [z],
                text  = string(node.ID),
                color = P[:node_label_color ],
                align = P[:align_node_labels])
        end
    end

    # Return the updated model plot:
    return P
end