function Base.show(io::IO, node::Node)
    print(io, styled"""
    {bold, cyan:Node:}
    ├─ x: $(node.x)
    ├─ y: $(node.y)
    └─ z: $(node.z)""")

    return nothing
end

function Base.show(io::IO, material::Material)
    print(io, styled"""
    {bold, cyan:Material:}
    ├─ E: $(material.E)
    ├─ ν: $(material.ν)
    └─ ρ: $(material.ρ)""")

    return nothing
end

function Base.show(io::IO, section::Section)
    print(io, styled"""
    {bold, cyan:Section:}
    ├─ A:    $(section.A   )
    ├─ I_zz: $(section.I_zz)
    ├─ I_yy: $(section.I_yy)
    └─ J:    $(section.J   )""")

    return nothing
end

function Base.show(io::IO, element::Element)
    print(io, styled"""
    {bold, cyan:Element ID:}
    ├─ {cyan:Node (i):}
    │  ├─ ID: $(element.node_i_ID)
    │  ├─ x:  $(element.x_i      )
    │  ├─ y:  $(element.y_i      )
    │  └─ z:  $(element.z_i      )
    ├─ {cyan:Node (j):}
    │  ├─ ID: $(element.node_j_ID)
    │  ├─ x:  $(element.x_j      )
    │  ├─ y:  $(element.y_j      )
    │  └─ z:  $(element.z_j      )
    ├─ {cyan:Material:}
    │  ├─ ID: $(element.material_ID)
    │  ├─ E:  $(element.E          )
    │  ├─ ν:  $(element.ν          )
    │  └─ ρ:  $(element.ρ          )
    ├─ {cyan:Section:}
    │  ├─ ID:   $(element.section_ID)
    │  ├─ A:    $(element.A         )
    │  ├─ I_zz: $(element.I_zz      )
    │  ├─ I_yy: $(element.I_yy      )
    │  └─ J:    $(element.J         )
    └─ ω: $(element.ω)""")

    return nothing
end

function Base.show(io::IO, support::Support)
    print(io, styled"""
    {bold, cyan:Support:}
    ├─ u_x constrained: $(support.u_x)
    ├─ u_y constrained: $(support.u_y)
    ├─ u_z constrained: $(support.u_z)
    ├─ θ_x constrained: $(support.θ_x)
    ├─ θ_y constrained: $(support.θ_y)
    └─ θ_z constrained: $(support.θ_z)""")

    return nothing
end

function Base.show(io::IO, cload::CLoad)
    print(io, styled"""
    {bold, cyan:Concentrated load:}
    ├─ F_x: $(cload.F_x)
    ├─ F_y: $(cload.F_y)
    ├─ F_z: $(cload.F_z)
    ├─ M_x: $(cload.M_x)
    ├─ M_y: $(cload.M_y)
    └─ M_z: $(cload.M_z)""")

    return nothing
end

function Base.show(io::IO, model::Model)
    num_nodes      = length(model.nodes    )
    num_sections   = length(model.sections )
    num_materials  = length(model.materials)
    num_elements   = length(model.elements )
    num_supports   = length(model.supports )
    num_cloads     = length(model.cloads   )
    num_dloads     = length(model.dloads   )

    if num_nodes == 0 && num_sections == 0 && num_materials == 0 && num_elements == 0 && num_supports == 0 && num_cloads == 0 && num_dloads == 0
        println(io, styled"""
        {cyan, bold:Empty model.}""")
    else
        println(io, styled"""
        {cyan, bold:Model:}
        ├─ # nodes      : $(num_nodes    )
        ├─ # sections   : $(num_materials)
        ├─ # materials  : $(num_sections )
        ├─ # elements   : $(num_elements )
        ├─ # supports   : $(num_supports )
        ├─ # conc. loads: $(num_cloads   )
        └─ # dist. loads: $(num_dloads   )""")
    end

    return nothing
end