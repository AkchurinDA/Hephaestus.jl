module Hephaestus
using OrderedCollections
using StyledStrings
using DocStringExtensions
using LinearAlgebra

include("Nodes.jl")
include("Materials.jl")
include("Sections.jl")
include("Elements.jl")
include("Models.jl")
include("ModelUtilities.jl")
include("Analyses/Analyses.jl")
export Model
export add_node!, add_material!, add_section!, add_element!, add_support!, add_nodal_load!
export del_node!, del_material!, del_section!, del_element!, del_support!, add_nodal_load!
export O1E, O2E
export solve
end