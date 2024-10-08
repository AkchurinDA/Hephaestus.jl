module HephaestusMakieExtension
using Hephaestus
using Makie

import Hephaestus: plotmodel, plotmodel!
import Hephaestus: plotsolution, plotsolution!

include("PlotModel.jl")
include("PlotSolution.jl")
end