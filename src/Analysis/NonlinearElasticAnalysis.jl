struct NonlinearElasticAnalysis <: AbstractAnalysisType
    
end

struct NonlinearElasticAnalysisCache{
    UT <: Real,
    RT <: Real} <: AbstractSolutionCache
    U::AbstractVector{UT}
    R::AbstractVector{RT}
    planar::Bool
end

function solve(model::Model, analysis::NonlinearElasticAnalysis, partitionindices::Vector{Bool})
    
end