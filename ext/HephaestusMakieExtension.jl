module HephaestusMakieExtension
using Makie
using Hephaestus

import Hephaestus: plotmodel, plotmodel!

Makie.@recipe(PlotModel, Model) do scene
    Makie.Attributes(
        # Nodes:
        nodevisible     = true,
        nodecolor       = :blue,
        nodemarker      = :circle,
        nodemarkersize  = 6,
        nodestrokecolor = :black,
        nodestrokewidth = 1,

        # Node labels:
        nodelabelvisible = true,
        nodelabelalign   = (:left, :bottom),
        nodelabelcolor   = :blue,

        # Elements:
        elementvisible   = true,
        elementcolor     = :black,
        elementlinestyle = :solid,
        elementlinewidth = 1,

        # Element labels:
        elementlabelvisible = true,
        elementlabelalign   = (:right, :top),
        elementlabelcolor   = :black,

        # Loads:
        concloadsvisible = true,
        concloadcolor    = :red,
        distloadsvisible = true,
        distloadcolor    = :red)
end

# Define the preferred axis type for the model plot:
function Makie.preferred_axis_type(P::PlotModel)
    if P[:Model][].dimensionality == 2
        return Makie.Axis
    end

    if P[:Model][].dimensionality == 3
        return Makie.Axis3
    end
end

# Define the plotting function for the model:
function Makie.plot!(P::PlotModel)
    # Extract the model:
    model = P[:Model]

    # Plot the elements:
    if P[:elementvisible][]
        if !isempty(model[].elements)
            for element in values(model[].elements)
                # Extract the coordinates of the element's nodes:
                x_i, y_i, z_i = element.node_i.x, element.node_i.y, element.node_i.z
                x_j, y_j, z_j = element.node_j.x, element.node_j.y, element.node_j.z

                # Plot the element:
                lines!(P, [x_i, x_j], [y_i, y_j], [z_i, z_j],
                    color     = P[:elementcolor    ],
                    linestyle = P[:elementlinestyle],
                    linewidth = P[:elementlinewidth])

                # Plot the element label:
                if P[:elementlabelvisible][]
                    # Compute the midpoint of the element:
                    x_m = (x_i + x_j) / 2
                    y_m = (y_i + y_j) / 2
                    z_m = (z_i + z_j) / 2

                    text!(P, x_m, y_m, z_m,
                        text  = "E$(element.ID)",
                        color = P[:elementlabelcolor],
                        align = P[:elementlabelalign])
                end
            end
        end
    else
        @warn "The model has no elements to plot."
    end

    # Plot the nodes:
    if P[:nodevisible][]
        if !isempty(model[].nodes)
            for node in values(model[].nodes)
                # Extract the coordinates of the node:
                x, y, z = node.x, node.y, node.z

                # Plot the node:
                scatter!(P, x, y, z,
                    color       = P[:nodecolor      ],
                    marker      = P[:nodemarker     ],
                    markersize  = P[:nodemarkersize ],
                    strokecolor = P[:nodestrokecolor],
                    strokewidth = P[:nodestrokewidth])

                # Plot the node label:
                if P[:nodelabelvisible][]
                    text!(P, [x], [y], [z],
                        text  = "N$(node.ID)",
                        color = P[:nodelabelcolor],
                        align = P[:nodelabelalign])
                end
            end
        end
    else
        @warn "The model has no nodes to plot."
    end

    # Plot the concentrated loads:
    if P[:concloadsvisible][]
        if !isempty(model[].concloads)

        end
    else
        @warn "The model has no concentrated loads to plot."
    end

    # Plot the distributed loads:
    if P[:distloadsvisible][]
        if !isempty(model[].distloads)

        end
    else
        @warn "The model has no distributed loads to plot."
    end

    # Return the updated model plot:
    return P
end
end
