abstract type AbstractModel end

@kwdef mutable struct Model
    materials = []
    sections  = []
    nodes     = []
    elements  = []
end
