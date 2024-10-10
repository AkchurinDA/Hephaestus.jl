"""
    struct O2EAnalysis

A type that represents the 2nd-order (`O2`) elastic (`E`) analysis.
"""
struct O2EAnalysis <: AbstractAnalysisType
    steps        ::Int
end

"""
    struct O2ESolutionCache

A type that stores the results of the 2nd-order elastic analysis.
"""
struct O2ESolutionCache <: AbstractSolutionCache

end


function solve(model::Model, analysis::O2EAnalysis)::O2ESolutionCache
    
end