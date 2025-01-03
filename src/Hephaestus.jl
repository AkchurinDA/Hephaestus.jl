module Hephaestus
using LinearAlgebra
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
export LinearElasticAnalysis, LinearElasticAnalysisCache
export NonlinearElasticAnalysis, NonlinearElasticAnalysisCache
export ElasticBucklingAnalysis, ElasticBucklingAnalysisCache
export FreeVibrationAnalysis, FreeVibrationAnalysisCache
export LCM, DCM, ALCM, WCM
export solve, solve!
export getnodaldisplacements, getnodalreactions
export getelementdisplacements, getelementforces
export generatereport
export plotmodel, plotmodel!
end
