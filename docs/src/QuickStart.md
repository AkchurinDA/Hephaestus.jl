To quickly get you started with `Hephaestus.jl`, let us recreate a simple example of a cantilever beam subjected to a concentrated load at its free end as shown below.

```@raw html

```

```@setup QuickStart
using Hephaestus
```

To create a new model, use the [`Model()`](@ref) constructor:

```@example QuickStart
M = Model()
```

To add nodes to the model, use the [`add_node!()`](@ref) function:

```@example QuickStart
add_node!(M, 1 , 0  * 12, 0, 0)
add_node!(M, 2 , 1  * 12, 0, 0)
add_node!(M, 3 , 2  * 12, 0, 0)
add_node!(M, 4 , 3  * 12, 0, 0)
add_node!(M, 5 , 4  * 12, 0, 0)
add_node!(M, 6 , 5  * 12, 0, 0)
add_node!(M, 7 , 6  * 12, 0, 0)
add_node!(M, 8 , 7  * 12, 0, 0)
add_node!(M, 9 , 8  * 12, 0, 0)
add_node!(M, 10, 9  * 12, 0, 0)
add_node!(M, 11, 10 * 12, 0, 0)
```

To add materials to the model, use the [`add_material!()`](@ref) function:

```@example QuickStart
add_material!(M, 1, 29000, 0.3, 0.000290)
```

To add sections to the model, use the [`add_section!()`](@ref) function:

```@example QuickStart
add_section!(M, 1, 9.16, 180, 180, 359)
```

To add elements to the model, use the [`add_element!()`](@ref) function:

```@example QuickStart
add_element!(M, 1 , 1 , 2 , 1, 1)
add_element!(M, 2 , 2 , 3 , 1, 1)
add_element!(M, 3 , 3 , 4 , 1, 1)
add_element!(M, 4 , 4 , 5 , 1, 1)
add_element!(M, 5 , 5 , 6 , 1, 1)
add_element!(M, 6 , 6 , 7 , 1, 1)
add_element!(M, 7 , 7 , 8 , 1, 1)
add_element!(M, 8 , 8 , 9 , 1, 1)
add_element!(M, 9 , 9 , 10, 1, 1)
add_element!(M, 10, 10, 11, 1, 1)
```

