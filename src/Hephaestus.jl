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
export Model
export Node, addnode!
export Material, addmaterial!
export Section, addsection!
export Element, addelement!
end
