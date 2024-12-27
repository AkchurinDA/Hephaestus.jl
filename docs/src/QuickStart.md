# Quick Start

## Defining a model

To quickly get you started with `Hephaestus.jl`, let us recreate a simple example of a cantilever beam subjected to a distributed load as shown below.

First of all, load the package using the following command:

```@example QuickStart
using Hephaestus
```

To create a new model, use the [`Model()`](@ref) constructor:

```@example QuickStart
model = Model()
```

To add nodes to the model, use the [`node!()`](@ref) function:

```@example QuickStart
node!(model,  1,  0 * 12, 0, 0, u_x = true, u_y = true, θ_z = true)
node!(model,  2,  1 * 12, 0, 0)
node!(model,  3,  2 * 12, 0, 0)
node!(model,  4,  3 * 12, 0, 0)
node!(model,  5,  4 * 12, 0, 0)
node!(model,  6,  5 * 12, 0, 0)
node!(model,  7,  6 * 12, 0, 0)
node!(model,  8,  7 * 12, 0, 0)
node!(model,  9,  8 * 12, 0, 0)
node!(model, 10,  9 * 12, 0, 0)
node!(model, 11, 10 * 12, 0, 0)
```

To add sections to the model, use the [`section!()`](@ref) function:

```@example QuickStart
section!(model, 1, 9.16, 180, 180, 359)
```

To add materials to the model, use the [`material!()`](@ref) function:

```@example QuickStart
material!(model, 1, 29000, 0.3, 0.290)
```

To add elements to the model, use the [`element!()`](@ref) function:

```@example QuickStart
element!(model, 1 , 1 , 2 , 1, 1)
element!(model, 2 , 2 , 3 , 1, 1)
element!(model, 3 , 3 , 4 , 1, 1)
element!(model, 4 , 4 , 5 , 1, 1)
element!(model, 5 , 5 , 6 , 1, 1)
element!(model, 6 , 6 , 7 , 1, 1)
element!(model, 7 , 7 , 8 , 1, 1)
element!(model, 8 , 8 , 9 , 1, 1)
element!(model, 9 , 9 , 10, 1, 1)
element!(model, 10, 10, 11, 1, 1)
```

To add distributed loads to the model, use the [`distload!()`](@ref) function:

```@example QuickStart
distload!(model,  1, 0, -1, 0)
distload!(model,  2, 0, -1, 0)
distload!(model,  3, 0, -1, 0)
distload!(model,  4, 0, -1, 0)
distload!(model,  5, 0, -1, 0)
distload!(model,  6, 0, -1, 0)
distload!(model,  7, 0, -1, 0)
distload!(model,  8, 0, -1, 0)
distload!(model,  9, 0, -1, 0)
distload!(model, 10, 0, -1, 0)
```

## Performing analysis

## Plotting model and solution

