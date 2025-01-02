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
function solve!(model::Model, analysistype::AbstractAnalysisType; continueanalaysis::Bool = false)::Model
    # Extract the partition indices:
    partitionindices = getpartitionindices(model)

    if continueanalaysis ≠ true
        # Initialize the state of the model:
        @info "The state of the model has been initialized."

        # Return the node states to their initial values:
        for node in model.nodes
            initstate!(node)
        end

        # Return the element states to their initial values:
        for element in model.elements
            initstate!(element)
        end
    end

    # Solve the model using the specified analysis type:
    solve!(model, analysistype, partitionindices)

    # Return the updated model:
    return model
end
