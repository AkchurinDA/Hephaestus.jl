module Hephaestus
using OrderedCollections
using LinearAlgebra
using StyledStrings

using DocStringExtensions

include("Nodes.jl")
export AbstractNode, Node
include("Materials.jl")
export AbstractMaterial, Material
include("Sections.jl")
export AbstractSection, Section
include("Elements.jl")
export AbstractElement, Element
include("Models.jl")
export AbstractModel, Model
export add_node!, add_material!, add_section!, add_element!, add_support!, add_conc_load!, add_dist_load!
export get_node_u_g, get_element_u_l, get_element_f_l
include("Analysis.jl")
export AbstractAnalysisType, AbstractAnalysisCache
export LinearElasticAnalysis, LinearElasticAnalysisCache
export NonlinearElasticAnalysis, NonlinearElasticAnalysisCache
export ElasticBucklingAnalysis, ElasticBucklingAnalysisCache
export FreeVibrationAnalysis, FreeVibrationAnalysisCache
export solve
include("Plotting.jl")
export plotmodel, plotmodel!
export plotsolution, plotsolution!
end