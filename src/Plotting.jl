"""
    plotmodel(model::Model, [options])

Plots a model of a structure of interest ([`Model`](@ref)) into a new `Makie.jl` scene.

# Plotting Options

## Nodes

| Option            | Description               | Default   |
| :---------------- | :------------------------ | :-------- |
| `node_visible`    | Whether to show the nodes | `true`    |
| `node_color`      | Color of the nodes        | `:red`    |
| `node_marker`     | Marker of the nodes       | `:circle` |
| `node_markersize` | Size of the nodes         | `6`       |
| `node_strokecolor`| Stroke color of the nodes | `:black`  |
| `node_strokewidth`| Stroke width of the nodes | `1`       |

## Node Labels

| Option               | Description                  | Default              |
| :------------------- | :--------------------------- | :------------------- |
| `node_label_visible` | Whether to label the nodes   | `true`               |
| `align_node_labels`  | Alignment of the node labels | `(:center, :bottom)` |
| `node_label_color`   | Color of the node labels     | `:red`               |

## Elements

| Option              | Description                  | Default  |
| :------------------ | :--------------------------- | :------- |
| `element_visible`   | Whether to show the elements | `true`   |
| `element_color`     | Color of the elements        | `:black` |
| `element_linestyle` | Line style of the elements   | `:solid` |
| `element_linewidth` | Line width of the elements   | `1`      |

## Element Labels

| Option                  | Description                     | Default              |
| :---------------------- | :------------------------------ | :------------------- |
| `element_label_visible` | Whether to label the elements   | `true`               |
| `align_element_labels`  | Alignment of the element labels | `(:center, :bottom)` |
| `element_label_color`   | Color of the element labels     | `:black`             |

## Supports

| Option            | Description                  | Default  |
| :---------------- | :--------------------------- | :------- |
| `support_visible` | Whether to show the supports | `true`   |
| `support_color`   | Color of the supports        | `:green` |
"""
function plotmodel end

"""
    plotmodel!(model::Model, [options])

Plots a model of a structure of interest ([`Model`](@ref)) into an existing scene. 
The plotting options are the same as in [`plotmodel()`](@ref) function.
"""
function plotmodel! end

"""
    plotsolution(model::Model, solution::AbstractSolutionCache, [options])

Plots undisplaced and displaced shapes of the model of a structure of interest ([`Model`](@ref)) into a new `Makie.jl` scene.

# Plotting Options

| Option  | Description  | Default |
| :------ | :----------- | :------ |
| `scale` | Scale factor | `1`     |

## Nodes of the Undisplaced Shape

| Option                    | Description               | Default      |
| :------------------------ | :-------------------------| :----------- |
| `undisp_node_visible`     | Whether to show the nodes | `true`       |
| `undisp_node_color`       | Color of the nodes        | `:steelblue` |
| `undisp_node_marker`      | Marker of the nodes       | `:circle`    |
| `undisp_node_markersize`  | Size of the nodes         | `6`          |
| `undisp_node_strokecolor` | Stroke color of the nodes | `:black`     |
| `undisp_node_strokewidth` | Stroke width of the nodes | `1`          |

## Node Labels of the Undisplaced Shape

| Option                      | Description                  | Default              |
| :-------------------------- | :--------------------------- | :------------------- |
| `undisp_node_label_visible` | Whether to label the nodes   | `true`               |
| `undisp_node_label_align`   | Alignment of the node labels | `(:center, :bottom)` |
| `undisp_node_label_color`   | Color of the node labels     | `:steelblue`         |

## Elements of the Undisplaced Shape

| Option                     | Description                  | Default  |
| :------------------------- | :--------------------------- | :------- |
| `undisp_element_visible`   | Whether to show the elements | `true`   |
| `undisp_element_color`     | Color of the elements        | `:black` |
| `undisp_element_linestyle` | Line style of the elements   | `:dash`  |
| `undisp_element_linewidth` | Line width of the elements   | `1`      |

## Nodes of the Displaced Shape

| Option                  | Description               | Default   |
| :---------------------- | :------------------------ | :-------- |
| `disp_node_visible`     | Whether to show the nodes | `true`    |
| `disp_node_color`       | Color of the nodes        | `:crimson`|
| `disp_node_marker`      | Marker of the nodes       | `:circle` |
| `disp_node_markersize`  | Size of the nodes         | `6`       |
| `disp_node_strokecolor` | Stroke color of the nodes | `:black`  |
| `disp_node_strokewidth` | Stroke width of the nodes | `1`       |

## Node Labels of the Displaced Shape

| Option                    | Description                  | Default           |
| :------------------------ | :--------------------------- | :---------------- |
| `disp_node_label_visible` | Whether to label the nodes   | `true`            |
| `disp_node_label_align`   | Alignment of the node labels | `(:center, :top)` |
| `disp_node_label_color`   | Color of the node labels     | `:crimson`        |

## Elements of the Displaced Shape

| Option                   | Description                  | Default  |
| :----------------------- | :--------------------------- | :------- |
| `disp_element_visible`   | Whether to show the elements | `true`   |
| `disp_element_color`     | Color of the elements        | `:black` |
| `disp_element_linestyle` | Line style of the elements   | `:solid` |
| `disp_element_linewidth` | Line width of the elements   | `1`      |
"""
function plotsolution end

"""
    plotsolution!(model::Model, [options])

Plots undisplaced and displaced shapes of the model of a structure of interest ([`Model`](@ref)) into an existing scene. 
The plotting options are the same as in [`plotsolution()`](@ref) function.
"""
function plotsolution! end