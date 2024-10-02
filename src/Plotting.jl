"""
    plotmodel(model::Model, [options])

Plots a model of a structure of interest ([`Model`](@ref)) into a new `Makie.jl` scene.

# Plotting Options

## Nodes

| Option            | Description               | Default   |
| :---------------- | :------------------------ | :-------- |
| `show_nodes`      | Whether to show the nodes | `true`    |
| `node_color`      | Color of the nodes        | `:red`    |
| `node_marker`     | Marker of the nodes       | `:circle` |
| `node_markersize` | Size of the nodes         | `6`       |
| `node_strokecolor`| Stroke color of the nodes | `:black`  |
| `node_strokewidth`| Stroke width of the nodes | `1`       |

## Node Labels

| Option             | Description                  | Default              |
| :----------------- | :--------------------------- | :------------------- |
| `label_nodes`      | Whether to label the nodes   | `true`               |
| `align_node_labels`| Alignment of the node labels | `(:center, :bottom)` |
| `node_label_color` | Color of the node labels     | `:red`               |

## Elements

| Option              | Description                  | Default  |
| :------------------ | :--------------------------- | :------- |
| `show_elements`     | Whether to show the elements | `true`   |
| `element_color`     | Color of the elements        | `:black` |
| `element_linestyle` | Line style of the elements   | `:solid` |
| `element_linewidth` | Line width of the elements   | `1`      |

## Element Labels

| Option                | Description                     | Default              |
| :-------------------- | :------------------------------ | :------------------- |
| `label_elements`      | Whether to label the elements   | `true`               |
| `align_element_labels`| Alignment of the element labels | `(:center, :bottom)` |
| `element_label_color` | Color of the element labels     | `:black`             |

## Supports

| Option          | Description                  | Default  |
| :-------------- | :--------------------------- | :------- |
| `show_supports` | Whether to show the supports | `true`   |
| `support_color` | Color of the supports        | `:green` |
"""
function plotmodel  end

"""
    plotmodel!(model::Model, [options])

Plots a model of a structure of interest ([`Model`](@ref)) into an existing scene. The plotting options are the same as in [`plotmodel()`](@ref) function.
"""
function plotmodel! end