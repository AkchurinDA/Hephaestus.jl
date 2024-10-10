<img src="docs/src/assets/social-preview.svg" alt="Hephaestus.jl">

<div align = "center">

  | Developer | [Damir Akchurin](https://scholar.google.com/citations?user=chYaDcIAAAAJ&hl=en) |
  | :--- | :--- |
  | Latest Release | [![Laterst Release](https://juliahub.com/docs/General/Hephaestus/0.1.0/version.svg)](https://juliahub.com/ui/Packages/General/Hephaestus) |
  | Documentation | [![Documentation](https://img.shields.io/badge/Documentation-Stable-blue.svg)](https://AkchurinDA.github.io/Hephaestus.jl/stable) <br> [![Documentation](https://img.shields.io/badge/Documentation-Dev-blue.svg)](https://AkchurinDA.github.io/Hephaestus.jl/dev) |
  | Downloads | [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FHephaestus&query=total_requests&label=Total)](http://juliapkgstats.com/pkg/Hephaestus) <br> [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FHephaestus&query=total_requests&label=Monthly&suffix=%2FMonth)](http://juliapkgstats.com/pkg/Hephaestus) |
  | License | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/AkchurinDA/Hephaestus.jl/blob/main/LICENSE.md) |

</div>

## Description

`Hephaestus.jl` is an auto-differentiable structural analysis package purely written in the Julia programming language.

## Installation

To install `Hephaestus.jl` package, type `]` in Julia REPL to enter the built-in Julia package manager and execute the following command:

```
pkg> add Hephaestus
```

## License

`Hephaestus.jl` package is distributed under the [MIT license](https://en.wikipedia.org/wiki/MIT_License). More information can be found in the [`LICENSE.md`](https://github.com/AkchurinDA/Hephaestus.jl/blob/main/LICENSE.md) file.

## Help and Support

For assistance with the package, please raise an issue on the [GitHub Issues](https://github.com/AkchurinDA/Hephaestus.jl/issues) page. Please use the appropriate labels to indicate the specific functionality you are inquiring about. Alternatively, contact the author directly at [AkchurinDA@gmail.com](mailto:AkchurinDA@gmail.com?subject=Hephaestus.jl).

## Acknowledgements

The design of the package is inspired by [`OpenSeesPy`](https://github.com/zhuminjie/OpenSeesPy), [`PyNite`](https://github.com/JWock82/Pynite), and [`MASTAN2`](https://www.mastan2.com).

## Roadmap

- [ ] Analyses
  - [x] 1nd-order elastic analysis
  - [ ] 2nd-order elastic analysis
  - [x] Elastic buckling analysis
  - [x] Free vibration analysis
- [ ] Elements
  - [ ] Truss element
  - [x] Beam-column element (Euler-Bernoulli)
  - [ ] Beam-column element (Timoshenko)
- [ ] Utilities
  - [x] Plotting a model
  - [ ] Extracting element information (displacement and force vectors in the local coordinate systems) from the solution cache