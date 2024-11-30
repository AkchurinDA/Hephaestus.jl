module Hephaestus
using OrderedCollections
using StyledStrings

abstract type AbstractNode     end
abstract type AbstractMaterial end
abstract type AbstractSection  end
abstract type AbstractElement  end

include("Components/Node.jl"    )
include("Components/Material.jl")
include("Components/Section.jl" )
include("Components/Element.jl" )
include("Components/Support.jl" )
include("Components/Load.jl"    )

include("Model.jl")

export Model
export add_node!, add_material!, add_section!, add_element!, add_support!, add_cload!, add_dload!
export del_node!, del_material!, del_section!, del_element!, del_support!, del_cload!, del_dload!

abstract type AbstractAnalysisType end
abstract type AbstractAnalysis     end

include("Analysis/LinearElasticAnalysis.jl"   )
include("Analysis/NonlinearElasticAnalysis.jl")
include("Analysis/ElasticBucklingAnalysis.jl" )
include("Analysis/FreeVibrationAnalysis.jl"   )

include("Analysis/Common.jl")

export LinearElasticAnalysis   , LinearElasticAnalysisCache
export NonlinearElasticAnalysis, NonlinearElasticAnalysisCache
export ElasticBucklingAnalysis , ElasticBucklingAnalysisCache
export FreeVibrationAnalysis   , FreeVibrationAnalysisCache
export solve

include("Utilities/Printing.jl")
end
