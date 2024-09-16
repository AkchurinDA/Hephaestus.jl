module Hephaestus
using SparseArrays
using StyledStrings
using OrderedCollections
using DocStringExtensions

include("Node.jl")
include("Material.jl")
include("Section.jl")
include("Element.jl")
include("Model.jl")
include("Analysis.jl")
export Node, Material, Section, Element, Model
export add_node!, add_material!, add_section!, add_element!, add_support!, add_nodal_load!, add_nodal_disp!
export del_node!, del_material!, del_section!, del_element!, del_support!, del_nodal_load!, del_nodal_disp!
export reset_model!
export solve
end