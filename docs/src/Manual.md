# Manual

```@setup Manual
using Hephaestus
```

## Defining a model

First, let's initialize a new model using the [`Model()`](@ref) constructor:

```@example Manual
M = Model()
```

## Adding nodes to the model

To add nodes to the model, use the [`add_node!()`](@ref) function:

```@example Manual
# add_node!(Model, Node ID, 
#   x-coordinate, y-coordinate, z-coordinate)
add_node!(M, 1, 0  , 0, 0)
add_node!(M, 2, 120, 0, 0)
```

## Adding materials to the model

To add materials to the model, use the [`add_material!()`](@ref) function:

```@example Manual
# add_material!(Model, Material ID, 
#   Young's modulus, Poisson's ratio, Density)
add_material!(M, 1, 29000, 0.3, 0.000290)
```

## Adding sections to the model

To add sections to the model, use the [`add_section!()`](@ref) function:

```@example Manual
# add_section!(Model, Section ID, 
#   Cross-sectional area, Moment of inertia about local z-axis, Moment of inertia about local y-axis, Polar moment of inertia)
add_section!(M, 1, 11.7, 307, 44.1, 0.91)
```