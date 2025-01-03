# News

## Release V0.3.0

- Dimensionality of a model is now carried as a parameter of `Model()` type.

## Release V0.2.0

- A completely new internal design of the package.
- Added separate types for conc. and dist. loads.
- Added `NodeState()` and `ElementState()` types to track the current states of the model. These types were needed to allow for nonlinear solvers to work.
- Added an option to easily run planar analysis using `solve(model, analysistype, planar = true)`, i.e., automatically restrict $u_{z}$, $\theta_{x}$, and $\theta_{y}$ DOFs by overwriting a user's input.
- Fully automatically differentiable structural analysis in now fully functional.
- Better pretty-printing of a model's components.
- Added automatic report generation that can be invoked by using `generatereport()` function.
- Added a plotting extension for quick model prototyping. This plotting feature is enabled by having an extension for all `Makie.jl`'s backends.
- Added more unit tests.

## Release V0.1.0

- Initial release of `Hephaestus.jl`.