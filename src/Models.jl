@kwdef struct Model
    nodes           ::OrderedDict{Int, Node    } = OrderedDict{Int, Node    }()
    materials       ::OrderedDict{Int, Material} = OrderedDict{Int, Material}()
    elements        ::OrderedDict{Int, Element } = OrderedDict{Int, Element }()
    sections        ::OrderedDict{Int, Section } = OrderedDict{Int, Section }()
    supports        ::OrderedDict{Int, Vector  } = OrderedDict{Int, Vector  }()
    nodal_loads     ::OrderedDict{Int, Vector  } = OrderedDict{Int, Vector  }()
end

function Base.show(io::IO, model::Model)
    if isempty(model.nodes) && isempty(model.materials) && isempty(model.elements) && isempty(model.sections)
        # Print the model information:
        printstyled(io, "Empty model \n", color = :yellow, bold = true)
    else 
        # Get the number of nodes, materials, elements, and sections:
        num_nodes     = length(model.nodes    )
        num_materials = length(model.materials)
        num_elements  = length(model.elements )
        num_sections  = length(model.sections )

        # Get the maximum number of digits in the number of nodes, materials, elements, and sections:
        str_length = round(Int, log10(max(num_nodes, num_materials, num_elements, num_sections)), RoundUp)

        # Stringify the number of nodes, materials, elements, and sections:
        str_num_nodes     = lpad(string(num_nodes    ), 1 + str_length, " ")
        str_num_materials = lpad(string(num_materials), 1 + str_length, " ")
        str_num_elements  = lpad(string(num_elements ), 1 + str_length, " ")
        str_num_sections  = lpad(string(num_sections ), 1 + str_length, " ")

        # Extract the tags of the nodes, materials, sections, and elements:
        # tags_nodes     = collect(keys(model.nodes    ))
        # tags_materials = collect(keys(model.materials))
        # tags_sections  = collect(keys(model.sections ))
        # tags_elements  = collect(keys(model.elements ))

        # Print the model information:
        printstyled(io, "Model information \n", color = :cyan, bold = true)
        printstyled(io, "└───$(str_num_nodes    ) \t Nodes     \n", color = :grey) #; printstyled(io, """  └───ID: $(join(tags_nodes    , ", ")) \n""", color = :grey)
        printstyled(io, "└───$(str_num_materials) \t Materials \n", color = :grey) #; printstyled(io, """  └───ID: $(join(tags_materials, ", ")) \n""", color = :grey)
        printstyled(io, "└───$(str_num_sections ) \t Sections  \n", color = :grey) #; printstyled(io, """  └───ID: $(join(tags_sections , ", ")) \n""", color = :grey)
        printstyled(io, "└───$(str_num_elements ) \t Elements  \n", color = :grey) #; printstyled(io, """  └───ID: $(join(tags_elements , ", ")) \n""", color = :grey) 
    end
end