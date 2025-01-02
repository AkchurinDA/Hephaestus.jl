# API

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
LinearElasticAnalysis
LinearElasticAnalysisCache
ElasticBucklingAnalysis
ElasticBucklingAnalysisCache
FreeVibrationAnalysis
FreeVibrationAnalysisCache
```

## Functions Used to Define a Model

```@docs
node!
section!
material!
element!
concload!
distload!
```

## Functions Used to Perform Analyses of Different Types and Extract the Results

```@docs
solve!
getnodaldisplacements
getnodalreactions
getelementdisplacements
getelementforces
```

## Functions Used to Plot a Model and the Results of Analyses of Different Types

```@docs
plotmodel
plotmodel!
```