struct NonlinearElasticAnalysis <: AbstractAnalysisType
    nonlinearsolver::AbstractNonlinearSolver
end

struct NonlinearElasticAnalysisCache{
    UT <: Real,
    RT <: Real} <: AbstractSolutionCache
    U::AbstractVector{UT}
    R::AbstractVector{RT}
end

function solve(model::Model, analysis::NonlinearElasticAnalysis, partitionindices::Vector{Bool})
    # Extract the nonlinear solver:
    nonlinearsolver = analysis.nonlinearsolver

    # Initialize the increment counter:
    for i in 1:maxnumi
        # Compute the initial tangent stiffness matrix:
        if j == 1 || update == :standard
            
        end

        # 
        if j == 1
            # δu_r = [0, 0, 0, ..., 0]
        else
            # δu_r = inv(K_t) * R
        end

        # Compute the load parameter increment:
        # δλ
    end
end