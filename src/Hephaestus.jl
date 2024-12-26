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

export Node, Section, Material, Element, ConcentratedLoad, DistributedLoad, Model
export node!, section!, material!, element!, concload!, distload!
export LinearElasticAnalysis, NonlinearElasticAnalysis
export solve
export extract_node_disp
export extract_element_disp_l, extract_element_force_l
export extract_element_disp_g, extract_element_force_g
export generatereport
end