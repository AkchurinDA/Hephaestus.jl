# Quick Start

## Defining a Model

To quickly get you started with `Hephaestus.jl`, let us recreate a simple example of a cantilever beam subjected to a distributed load as shown below.

```@raw html

```

First of all, load the package using the following command:

```@example QuickStart
using Hephaestus
```

To create a new model, use the [`Model()`](@ref) constructor:

```@example QuickStart
model = Model()
```

To add nodes to the model, use the [`add_node!()`](@ref) function:

```@example QuickStart
add_node!(model, 1 , 0  * 12, 0, 0)
add_node!(model, 2 , 1  * 12, 0, 0)
add_node!(model, 3 , 2  * 12, 0, 0)
add_node!(model, 4 , 3  * 12, 0, 0)
add_node!(model, 5 , 4  * 12, 0, 0)
add_node!(model, 6 , 5  * 12, 0, 0)
add_node!(model, 7 , 6  * 12, 0, 0)
add_node!(model, 8 , 7  * 12, 0, 0)
add_node!(model, 9 , 8  * 12, 0, 0)
add_node!(model, 10, 9  * 12, 0, 0)
add_node!(model, 11, 10 * 12, 0, 0)
```

To add materials to the model, use the [`add_material!()`](@ref) function:

```@example QuickStart
add_material!(model, 1, 29000, 0.3, 0.000290)
```

To add sections to the model, use the [`add_section!()`](@ref) function:

```@example QuickStart
add_section!(model, 1, 9.16, 180, 180, 359)
```

To add elements to the model, use the [`add_element!()`](@ref) function:

```@example QuickStart
add_element!(model, 1 , 1 , 2 , 1, 1)
add_element!(model, 2 , 2 , 3 , 1, 1)
add_element!(model, 3 , 3 , 4 , 1, 1)
add_element!(model, 4 , 4 , 5 , 1, 1)
add_element!(model, 5 , 5 , 6 , 1, 1)
add_element!(model, 6 , 6 , 7 , 1, 1)
add_element!(model, 7 , 7 , 8 , 1, 1)
add_element!(model, 8 , 8 , 9 , 1, 1)
add_element!(model, 9 , 9 , 10, 1, 1)
add_element!(model, 10, 10, 11, 1, 1)
```

To add supports to the model and fix particular degrees of freedom, use the [`add_support!()`](@ref) function:

```@example QuickStart
add_support!(model, 1 , true , true, true , true, true , true)
add_support!(model, 2 , false, true, false, true, false, true)
add_support!(model, 3 , false, true, false, true, false, true)
add_support!(model, 4 , false, true, false, true, false, true)
add_support!(model, 5 , false, true, false, true, false, true)
add_support!(model, 6 , false, true, false, true, false, true)
add_support!(model, 7 , false, true, false, true, false, true)
add_support!(model, 8 , false, true, false, true, false, true)
add_support!(model, 9 , false, true, false, true, false, true)
add_support!(model, 10, false, true, false, true, false, true)
add_support!(model, 11, false, true, false, true, false, true)
```

To add distributed loads to the model, use the [`add_dist_load!()`](@ref) function:

```@example QuickStart
add_dist_load!(model, 1 , 0, 0, -1)
add_dist_load!(model, 2 , 0, 0, -1)
add_dist_load!(model, 3 , 0, 0, -1)
add_dist_load!(model, 4 , 0, 0, -1)
add_dist_load!(model, 5 , 0, 0, -1)
add_dist_load!(model, 6 , 0, 0, -1)
add_dist_load!(model, 7 , 0, 0, -1)
add_dist_load!(model, 8 , 0, 0, -1)
add_dist_load!(model, 9 , 0, 0, -1)
add_dist_load!(model, 10, 0, 0, -1)
```

## Performing Analysis

To perform the first-order elastic analysis, use the [`solve()`](@ref) function with the second argument being [`LinearElasticAnalysis()`](@ref), which will make `Hephaestus.jl` to perform the first-order elastic analysis. Alternatively, to perform the:
- Second-order elastic analysis use [`NonlinearElasticAnalysis()`](@ref)
- Elastic buckling analysis use [`ElasticBucklingAnalysis()`](@ref)
- Free vibration analysis use [`FreeVibrationAnalysis()`](@ref)

```@example QuickStart
solution = solve(model, LinearElasticAnalysis())
```

## Plotting the Model and Solution

`Hephaestus.jl` exports two plotting recipes:

| Function | Description |
| :--- | :--- |
| [`plotmodel!()`](@ref) | Plots the defined model |
| [`plotsolution!()`](@ref) | Plots the undisplaced and displaced shapes of the model based on the solution cache |

```julia
using CairoMakie
CairoMakie.activate!(type = "svg")
CairoMakie.set_theme!(theme_latexfonts())

begin
    # Define the figure:
    F = Figure(size = 72 .* (12, 6))
   
    # Define the axis for the defined model:   
    A = Axis3(F[1, 1],
        title       = "Model",
        xlabel      = L"$x$ (in.)", xlabelrotation = 0,
        ylabel      = L"$y$ (in.)", ylabelrotation = 0,
        zlabel      = L"$z$ (in.)", zlabelrotation = +π / 2,
        limits      = (nothing, nothing, -60, +60, -60, +60),
        aspect      = :data,
        azimuth     = -π / 4, elevation = +π / 12,
        protrusions = 60)

    # Plot the defined model:
    plotmodel!(A, model)

    # Define the axis for the results of the first-order elastic analysis:
    A = Axis3(F[1, 2],
        title       = "Solution",
        xlabel      = L"$x$ (in.)", xlabelrotation = 0,
        ylabel      = L"$y$ (in.)", ylabelrotation = 0,
        zlabel      = L"$z$ (in.)", zlabelrotation = +π / 2,
        limits      = (nothing, nothing, -60, +60, -60, +60),
        aspect      = :data,
        azimuth     = -π / 4, elevation = +π / 12,
        protrusions = 60)

    # Plot the results of the first-order elastic analysis:
    plotsolution!(A, model, solution,
        scale = 5)

    # Display the figure:
    display(F)
end
```

```@raw html
<img src="../assets/QuickStart.pdf" class="center" style="max-height: 432px; border-radius:5px;"/>
```

