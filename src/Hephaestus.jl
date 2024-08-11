module Hephaestus
# --------------------------------------------------
# PREAMBLE
# --------------------------------------------------
import StaticArrays
import LinearAlgebra

# --------------------------------------------------
# TYPES
# --------------------------------------------------
include("Elements/Truss.jl")
export TrussElement
include("Elements/BeamColumnElement.jl")
export BeamColumnElement
end
