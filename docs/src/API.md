## Types Used to Define a Model

```@docs
Model
Node
Material
Section
Element
```

## Types Used to Perform Analyses of Different Types and Store the Results

```@docs
O1EAnalysis
O1ESolutionCache
O2EAnalysis
O2ESolutionCache
EBAnalysis
EBSolutionCache
FVAnalysis
FVSolutionCache
```

## Functions Used to Define a Model

```@docs
add_node!
del_node!
add_material!
del_material!
add_section!
del_section!
add_element!
del_element!
add_support!
del_support!
add_conc_load!
del_conc_load!
add_dist_load!
del_dist_load!
```

## Functions Used to Perform Analyses of Different Types and Extract the Results

```@docs
solve
get_node_u_g
get_element_u_l
get_element_f_l
```

## Functions Used to Plot a Model and the Results of Analyses of Different Types

```@docs
plotmodel
plotmodel!
plotsolution
plotsolution!
```