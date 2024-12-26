struct NonlinearElasticAnalysis <: AbstractAnalysisType
    
end

struct NonlinearElasticAnalysisCache{
    UT <: Real,
    RT <: Real} <: AbstractAnalysisCache
    U::AbstractVector{UT}
    R::AbstractVector{RT}
end

function solve(model::Model, analysis::NonlinearElasticAnalysis, indices_f::Vector{Bool}, indices_s::Vector{Bool})
    
end