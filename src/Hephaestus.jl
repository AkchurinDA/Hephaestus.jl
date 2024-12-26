module Hephaestus
using Memoization
using Dates
using DocStringExtensions
using Printf
using StyledStrings

include("Components/Nodes.jl")
include("Components/Sections.jl")
include("Components/Materials.jl")
include("Components/Elements.jl")
include("Components/Loads.jl")
include("Models.jl")
include("Analysis.jl")
include("Utilities/PrettyPrint.jl")
include("Utilities/GenerateReport.jl")
include("Utilities/Plotting.jl")

export Node, Section, Material, Element, ConcentratedLoad, DistributedLoad, Model
export node!, section!, material!, element!, concload!, distload!
export LinearElasticAnalysis, NonlinearElasticAnalysis
export solve
export getnodedisp
export getelementdisp_l, getelementforces_l
export getelementdisp_g, getelementforces_g
export generatereport
export plotmodel, plotmodel!
end