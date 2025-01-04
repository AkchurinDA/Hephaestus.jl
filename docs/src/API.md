# API

## Types

```@docs
Model
Node
Material
Section
Element
LinearElasticAnalysis
NonlinearElasticAnalysis
ElasticBucklingAnalysis
ElasticBucklingAnalysisCache
FreeVibrationAnalysis
FreeVibrationAnalysisCache
LCM
DCM
WCM
ALCM
```

## Constructing a model

```@docs
node!
section!
material!
element!
concload!
distload!
```

## Analyzing a model

```@docs
solve!
solve
```

## Postprocessing results

```@docs
getnodaldisplacements
getnodalreactions
getelementdisplacements
getelementforces
```

## Plotting a model

```@docs
plotmodel
plotmodel!
```