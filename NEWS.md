# News

## Release V0.2.0

- A completely new internal design of the package.
- Added separate types for conc. and dist. loads.
- Added an option to easily run planar analysis using `solve(model, analysistype, planar = true)`, i.e., automatically restrict $u_{z}$, $\theta_{x}$, and $\theta_{y}$ DOFs.
- Added a plotting extension for quick model prototyping.
- Better pretty-printing of a model's components.
- Added automatic report generation that can be invoked by using `generatereport()` function.

## Release V0.1.0

- Initial release of `Hephaestus.jl`.