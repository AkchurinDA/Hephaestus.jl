module Hephaestus
# --------------------------------------------------
# PREAMBLE
# --------------------------------------------------
import StaticArrays
import SparseArrays
import LinearAlgebra

# --------------------------------------------------
# TYPES
# --------------------------------------------------
include("Model.jl")
export 
    Model, 
    Node,
    Material, addmaterial!
    Section, addsection!
    Element, addelement!
    addnode!,
    addmaterial!,
    addsection!,
    addelement!
end
