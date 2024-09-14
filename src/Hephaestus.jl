module Hephaestus
    using SparseArrays
    
    include("Model.jl")
    export Node, Material, Section, Element, Model
    export add_node!, delete_node!, add_material!, delete_material!, add_section!, delete_section!, add_element!, delete_element!
    export clear_model!
end