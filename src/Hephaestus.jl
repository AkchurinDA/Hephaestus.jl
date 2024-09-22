module Hephaestus
using LinearAlgebra
using StyledStrings
using DocStringExtensions

include("Nodes.jl")
include("Materials.jl")
include("Sections.jl")
include("Elements.jl")
include("Models.jl")
include("Analyses.jl")
export Node, Material, Section, Element, Model
export add_node!, add_material!, add_section!, add_element!, add_support!, add_concentrated_load!, add_distributed_load!
export O1EAnalysis, O1ECache
export O2EAnalysis, O2ECache
export solve
end