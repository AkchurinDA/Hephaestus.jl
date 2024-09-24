module Hephaestus
using BlockDiagonals
using StyledStrings

include("Nodes.jl"    )
include("Materials.jl")
include("Sections.jl" )
include("Elements.jl" )
include("Models.jl"   )
include("Utilities.jl")
include("Analyses/Analyses.jl")
export Node, Material, Section, Element, Model
export O1EAnalysis, O1ECache
export add_node!, add_material!, add_section!, add_element!, add_support!, add_conc_load!, add_dist_load!
export del_node!, del_material!, del_section!, del_element!, del_support!, del_conc_load!, del_dist_load!
export solve
end