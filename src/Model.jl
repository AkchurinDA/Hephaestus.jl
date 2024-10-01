include("Node.jl")
include("Material.jl")
include("Section.jl")
include("Element.jl")
include("Spring.jl")

"""
    struct Model

A type representing the finite element model of a structure of interest.

$(FIELDS)
"""
@kwdef struct Model
    "Ordered dictionary that stores the nodes of the model"
    nodes         ::OrderedDict{Int, Node    } = OrderedDict{Int, Node    }()
    "Ordered dictionary that stores the materials of the model"
    materials     ::OrderedDict{Int, Material} = OrderedDict{Int, Material}()
    "Ordered dictionary that stores the sections of the model"
    sections      ::OrderedDict{Int, Section } = OrderedDict{Int, Section }()
    "Ordered dictionary that stores the elements of the model"
    elements      ::OrderedDict{Int, Element } = OrderedDict{Int, Element }()
    "Ordered dictionary that stores the springs of the model"
    springs       ::OrderedDict{Int, Spring  } = OrderedDict{Int, Spring  }()
    "Ordered dictionary that stores the supports of the model"
    supports      ::OrderedDict{Int, Vector{Bool}} = OrderedDict{Int, Vector{Bool}}()
    "Ordered dictionary that stores the concentrated loads of the model"
    conc_loads    ::OrderedDict{Int, Vector{<:Real}} = OrderedDict{Int, Vector{<:Real}}()
    "Ordered dictionary that stores the distributed loads of the model"
    dist_loads    ::OrderedDict{Int, Vector{<:Real}} = OrderedDict{Int, Vector{<:Real}}()
    "Ordered dictionary that stores the element fixed-end forces in the local coordinate system"
    p_l           ::OrderedDict{Int, Vector{<:Real}} = OrderedDict{Int, Vector{<:Real}}()
    "Ordered dictionary that stores the element fixed-end forces in the global coordinate system"
    p_g           ::OrderedDict{Int, Vector{<:Real}} = OrderedDict{Int, Vector{<:Real}}()
end

function Base.show(io::IO, model::Model)
    num_nodes      = length(model.nodes     )
    num_materials  = length(model.materials )
    num_sections   = length(model.sections  )
    num_elements   = length(model.elements  )
    num_supports   = length(model.supports  )
    num_conc_loads = length(model.conc_loads)
    num_dist_loads = length(model.dist_loads)

    if num_nodes == 0 && num_materials == 0 && num_sections == 0 && num_elements == 0 && num_supports == 0 && num_conc_loads == 0 && num_dist_loads == 0
        return println(io, styled"{cyan, bold: Empty model.}")
    else
        reserved_spacing = floor(Int, log10(maximum([num_nodes, num_materials, num_sections, num_elements, num_supports, num_conc_loads, num_dist_loads])))  + 1

        str_num_nodes      = lpad(num_nodes     , reserved_spacing)
        str_num_materials  = lpad(num_materials , reserved_spacing)
        str_num_sections   = lpad(num_sections  , reserved_spacing)
        str_num_elements   = lpad(num_elements  , reserved_spacing)
        str_num_supports   = lpad(num_supports  , reserved_spacing)
        str_num_conc_loads = lpad(num_conc_loads, reserved_spacing)
        str_num_dist_loads = lpad(num_dist_loads, reserved_spacing)

        println(io, styled"{cyan, bold: Model with:}")
        num_nodes      == 0 ? nothing : println(io, styled"{cyan: $(str_num_nodes     ) \t Nodes          }")
        num_materials  == 0 ? nothing : println(io, styled"{cyan: $(str_num_materials ) \t Materials      }")
        num_sections   == 0 ? nothing : println(io, styled"{cyan: $(str_num_sections  ) \t Sections       }")
        num_elements   == 0 ? nothing : println(io, styled"{cyan: $(str_num_elements  ) \t Elements       }")
        num_supports   == 0 ? nothing : println(io, styled"{cyan: $(str_num_supports  ) \t Supported nodes}")
        num_conc_loads == 0 ? nothing : println(io, styled"{cyan: $(str_num_conc_loads) \t Loaded nodes   }")
        num_dist_loads == 0 ? nothing : println(io, styled"{cyan: $(str_num_dist_loads) \t Loaded elements}")
    end

    return nothing
end