To quickly get you started with `Hephaestus.jl`, let us recreate a simple example of a cantilever beam subjected to a distributed load as shown below.

```@raw html

```

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

```@example QuickStart
solution = solve(model, O1EAnalysis())
```

```@example QuickStart
using CairoMakie
CairoMakie.activate!(type = "svg")
CairoMakie.set_theme!(theme_latexfonts()) # hide

begin
    F = Figure(size = 72 .* (8, 6))

    A = Axis3(F[1, 1],
        xlabel      = L"$x$ (in.)", xlabelrotation = 0,
        ylabel      = L"$y$ (in.)", ylabelrotation = 0,
        zlabel      = L"$z$ (in.)", zlabelrotation = +π / 2,
        limits      = (nothing, nothing, -60, +60, -60, +60),
        aspect      = :data,
        azimuth     = -π / 6, elevation = +π / 12)

    plotsolution!(A, model, solution,
        scale = 5)

    display(F)
end

save("QuickStart.pdf", F) # hide

nothing # hide
```

```@raw html
<img src="../QuickStart.pdf" class="center" style="max-width: 100%; border-radius:5px;"/>
```

