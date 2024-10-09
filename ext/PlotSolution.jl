Makie.@recipe(PlotSolution, Model, Solution) do scene
    Makie.Attributes(
        scale             = 1,

        undisp_node_visible     = true,
        undisp_node_color       = :steelblue,
        undisp_node_marker      = :circle,
        undisp_node_markersize  = 6,
        undisp_node_strokecolor = :black,
        undisp_node_strokewidth = 1,

        undisp_node_label_visible = false,
        undisp_node_label_align   = (:center, :bottom),
        undisp_node_label_color   = :steelblue,

        undisp_element_visible   = true,
        undisp_element_color     = :black,
        undisp_element_linestyle = :dash,
        undisp_element_linewidth = 1,

        disp_node_visible     = true,
        disp_node_color       = :crimson,
        disp_node_marker      = :circle,
        disp_node_markersize  = 6,
        disp_node_strokecolor = :black,
        disp_node_strokewidth = 1,

        disp_node_label_visible = false,
        disp_node_label_align   = (:center, :top),
        disp_node_label_color   = :crimson,

        disp_element_visible   = true,
        disp_element_color     = :black,
        disp_element_linestyle = :solid,
        disp_element_linewidth = 1)
end

Makie.preferred_axis_type(::PlotSolution) = Makie.Axis3

const PlotO1ESolution = PlotSolution{Tuple{<:Model, <:O1ESolutionCache}}
function Makie.plot!(P::PlotO1ESolution)
    model    = P[:Model]
    solution = P[:Solution]

    if !isempty(model[].elements)
        for element in values(model[].elements)
            x_i_old, y_i_old, z_i_old = element.x_i, element.y_i, element.z_i
            x_j_old, y_j_old, z_j_old = element.x_j, element.y_j, element.z_j 

            if P[:undisp_element_visible][]
                lines!(P, [x_i_old, x_j_old], [y_i_old, y_j_old], [z_i_old, z_j_old],
                    color     = P[:undisp_element_color    ],
                    linestyle = P[:undisp_element_linestyle],
                    linewidth = P[:undisp_element_linewidth])
            end

            if P[:disp_element_visible][]
                u_x_i, u_y_i, u_z_i, _, _, _ = get_node_u_g(solution[], element.node_i_ID)
                u_x_j, u_y_j, u_z_j, _, _, _ = get_node_u_g(solution[], element.node_j_ID)

                x_i_new = x_i_old + P[:scale][] * u_x_i
                y_i_new = y_i_old + P[:scale][] * u_y_i
                z_i_new = z_i_old + P[:scale][] * u_z_i
                x_j_new = x_j_old + P[:scale][] * u_x_j
                y_j_new = y_j_old + P[:scale][] * u_y_j
                z_j_new = z_j_old + P[:scale][] * u_z_j

                lines!(P, [x_i_new, x_j_new], [y_i_new, y_j_new], [z_i_new, z_j_new], 
                    color     = P[:disp_element_color    ],
                    linestyle = P[:disp_element_linestyle],
                    linewidth = P[:disp_element_linewidth])
            end
        end
    else
        @warn "The model has no elements to plot."
    end

    if !isempty(model[].nodes)
        for node in values(model[].nodes)
            x_old, y_old, z_old = node.x, node.y, node.z

            if P[:undisp_node_visible][]
                scatter!(P, [x_old], [y_old], [z_old],
                    color       = P[:undisp_node_color      ],
                    marker      = P[:undisp_node_marker     ],
                    markersize  = P[:undisp_node_markersize ],
                    strokecolor = P[:undisp_node_strokecolor],
                    strokewidth = P[:undisp_node_strokewidth],
                    overdraw    = true)

                if P[:undisp_node_label_visible][]
                    text!(P, [x_old], [y_old], [z_old],
                        text  = string(node.ID),
                        color = P[:undisp_node_label_color],
                        align = P[:undisp_node_label_align])
                end
            end

            if P[:disp_node_visible][]
                u_x, u_y, u_z, _, _, _ = get_node_u_g(solution[], node.ID)

                x_new = node.x + P[:scale][] * u_x
                y_new = node.y + P[:scale][] * u_y
                z_new = node.z + P[:scale][] * u_z

                scatter!(P, [x_new], [y_new], [z_new], 
                    color       = P[:disp_node_color      ],
                    marker      = P[:disp_node_marker     ],
                    markersize  = P[:disp_node_markersize ],
                    strokecolor = P[:disp_node_strokecolor],
                    strokewidth = P[:disp_node_strokewidth],
                    overdraw    = true)

                if P[:disp_node_label_visible][]
                    text!(P, [x_new], [y_new], [z_new],
                        text  = string(node.ID),
                        color = P[:disp_node_label_color],
                        align = P[:disp_node_label_align])
                end
            end
        end
    else
        @warn "The model has no nodes to plot."
    end
end