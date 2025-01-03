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
    solve!(model::Model, analysistype::AbstractAnalysisType)::Model

Function used to perform (geometrically) linear and nonlinear and (materially) elastic and inelastic analyses.

!!! note "Note"
    This function is mutating and will update the state of the model.
"""
function solve!(model::Model, analysistype::AbstractAnalysisType)::Model
    # This mutating method is only valid for geometrically linear and nonlinear and materially elastic and inelastic analyses:
    if !(isa(analysistype, LinearElasticAnalysis) || isa(analysistype, NonlinearElasticAnalysis))
        error("""
        The mutating `solve!()` method is only valid for (geometrically) linear and nonlinear and (materially) elastic analyses, i.e. `LinearElasticAnalysis()` and `NonlinearElasticAnalysis()`.
        To perform elastic buckling or free vibration analyses, use the nonmutating `solve()` method.""")
    end

    # Reinitialize the node states:
    for node in model.nodes
        initstate!(node)
    end

    # Reinitialize the element states:
    for element in model.elements
        initstate!(element)
    end

    # Extract the partition indices:
    partitionindices = getpartitionindices(model)

    # Solve the model using the specified analysis type:
    solve!(model, analysistype, partitionindices)

    # Return the updated model:
    return model
end

"""
    solve(model::Model, analysistype::AbstractAnalysisType)::AbstractSolutionCache

Function used to perform elastic buckling and free vibration analyses.

!!! note "Note"
    This function is nonmutating and will not update the state of the model. Instead, it return a solution cache object that can be used to extract the results.
"""
function solve(model::Model, analysistype::AbstractAnalysisType)::AbstractSolutionCache
    # This nonmutating method is only valid for elastic buckling and free vibration analyses:
    if !(isa(analysistype, ElasticBucklingAnalysis) || isa(analysistype, FreeVibrationAnalysis))
        error("""
        The nonmutating `solve()` method is only valid for elastic buckling or free vibration analyses, i.e. `ElasticBucklingAnalysis()` and `FreeVibrationAnalysis()`.
        To perform (geometrically) linear and nonlinear and (materially) elastic analyses, use the mutating `solve!()` method.""")
    end

    # Extract the partition indices:
    partitionindices = getpartitionindices(model)

    # Solve the model using the specified analysis type:
    solve(model, analysistype, partitionindices)

    # Return the updated model:
    return model
end
