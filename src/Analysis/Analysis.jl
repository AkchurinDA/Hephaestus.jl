"""
    solve(model::Model, analysis::AbstractAnalysisType)

The main function that performs the analysis of your choice on the model.
"""
function solve end

include("AnalysisUtilities.jl")
include("O1EAnalysis.jl")
include("O2EAnalysis.jl")
include("EBAnalysis.jl")
include("FVAnalysis.jl")