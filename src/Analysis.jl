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
function solve(model::Model, analysistype::AbstractAnalysisType)::AbstractSolutionCache
    # Extract the partition indices:
    partitionindices = getpartitionindices(model)

    # Return the node states to their initial values:
    if any([nodestate.modified for nodestate in model.nodestates]) # If any node state has been modified
        empty!(model.nodestates)
        for node in model.nodes
            push!(model.nodestates, initnodestate(node))
        end
        @info "Node states have been reset to their initial values."
    end

    # Return the element states to their initial values:
    if any([elementstate.modified for elementstate in model.elementstates]) # If any element state has been modified
        empty!(model.elementstates)
        for element in model.elements
            push!(model.elementstates, initelementstate(element))
        end
        @info "Element states have been reset to their initial values."
    end

    # Solve the model using the specified analysis type:
    solution = solve(model, analysistype, partitionindices)

    # Return the solution:
    return solution
end
