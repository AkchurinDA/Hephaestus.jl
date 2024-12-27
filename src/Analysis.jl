abstract type AbstractAnalysisType  end
abstract type AbstractSolutionCache end

include("Analysis/Assembling.jl")
include("Analysis/NonlinearSolvers.jl")
include("Analysis/LinearElasticAnalysis.jl")
include("Analysis/NonlinearElasticAnalysis.jl")
include("Analysis/ElasticBucklingAnalysis.jl")
include("Analysis/FreeVibrationAnalysis.jl")
include("Analysis/Postprocessing.jl")

"""
    solve(model::Model, analysistype::AbstractAnalysisType)

Solve the model using the specified analysis type.
"""
function solve(model::Model, analysistype::AbstractAnalysisType)
    # Extract the partition indices:
    partitionindices = getpartitionindices(model)

    # Solve the model using the specified analysis type:
    solution = solve(model, analysistype, partitionindices)

    # Return the solution:
    return solution
end
