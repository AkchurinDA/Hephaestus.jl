module HephaestusMakieExtension
using Hephaestus
using Makie

import Hephaestus: plotmodel, plotmodel!

Makie.@recipe(PlotModel, Model) do scene
    Makie.Attributes(
        # Nodes:
        show_nodes       = true,
        node_color       = :red ,
        node_marker      = :circle,
        node_markersize  = 6,
        node_strokecolor = :black,
        node_strokewidth = 1,

        # Node labels:
        show_node_labels = true,
        node_label_align = (:center, :bottom),
        node_label_color = :red,

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

    if P[:show_elements][] && !isempty(model[].elements)
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

    if P[:show_nodes][] && !isempty(model[].nodes)
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

    # if P[:show_supports][] && !isempty(model[].supports)
    #     largest_dimension_x = maximum([node.x for node in values(model[].nodes)])
    #     largest_dimension_y = maximum([node.y for node in values(model[].nodes)])
    #     largest_dimension_z = maximum([node.z for node in values(model[].nodes)])
    #     largest_dimension   = maximum([largest_dimension_x, largest_dimension_y, largest_dimension_z])
    #     support_size        = largest_dimension / 10

    #     for (ID, support) in model[].supports
    #         x, y, z = model[].nodes[ID].x, model[].nodes[ID].y, model[].nodes[ID].z

    #         u_x, u_y, u_z = support[1], support[2], support[3]

    #         if u_x
    #             arrows!(P, [x], [y], [z],
    #                 [support_size], [0], [0],
    #                 color     = P[:support_color],
    #                 linewidth = 0.1,
    #                 arrowsize = Vec3f(0.25, 0.25, 0.5),
    #                 shading   = NoShading)
    #         end

    #         if u_y
    #             arrows!(P, [x], [y], [z],
    #                 [0], [support_size], [0],
    #                 color     = P[:support_color],
    #                 linewidth = 0.1,
    #                 arrowsize = Vec3f(0.25, 0.25, 0.5),
    #                 shading   = NoShading)
    #         end

    #         if u_z
    #             arrows!(P, [x], [y], [z],
    #                 [0], [0], [support_size],
    #                 color     = P[:support_color],
    #                 linewidth = 0.1,
    #                 arrowsize = Vec3f(0.25, 0.25, 0.5),
    #                 shading   = NoShading)
    #         end
    #     end
    # end

    # Return the updated model plot:
    return P
end
end