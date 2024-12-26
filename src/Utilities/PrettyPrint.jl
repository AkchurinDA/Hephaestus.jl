Base.show(io::IO, node::Node) = print(io,
    styled"{cyan, bold:Node} with ID $(node.ID)")
Base.show(io::IO, ::MIME"text/plain", node::Node) = print(io, styled"""
    {cyan, bold:Node:}
    ├─ ID: $(node.ID)
    ├─ Nodal coordinates:
    │  ├─ x: $(node.x)
    │  ├─ y: $(node.y)
    │  └─ z: $(node.z)
    ├─ Translational DOF supports:
    │  ├─ $(node.u_x)
    │  ├─ $(node.u_y)
    │  └─ $(node.u_z)
    └─ Rotational DOF supports:
       ├─ $(node.θ_x)
       ├─ $(node.θ_y)
       └─ $(node.θ_z)""")

Base.show(io::IO, section::Section) = print(io,
    styled"{cyan, bold:Section} with ID $(section.ID)")
Base.show(io::IO, ::MIME"text/plain", section::Section) = print(io, styled"""
    {cyan, bold:Section:}
    ├─ ID:   $(section.ID  )
    ├─ A:    $(section.A   )
    ├─ I_zz: $(section.I_zz)
    ├─ I_yy: $(section.I_yy)
    └─ J:    $(section.J   )""")

Base.show(io::IO, material::Material) = print(io,
    styled"{cyan, bold:Material} with ID $(material.ID)")
Base.show(io::IO, ::MIME"text/plain", material::Material) = print(io, styled"""
    {cyan, bold:Material:}
    ├─ ID: $(material.ID)
    ├─ E:  $(material.E )
    ├─ ν:  $(material.ν )
    └─ ρ:  $(material.ρ )""")

Base.show(io::IO, element::Element) = print(io,
    styled"{cyan, bold:Element} with ID $(element.ID)")
Base.show(io::IO, ::MIME"text/plain", element::Element) = print(io, styled"""
    {cyan, bold:Element:}
    ├─ ID: $(element.ID)
    ├─ {bright_cyan:Node (i):}
    │  ├─ ID: $(element.node_i.ID)
    │  ├─ x:  $(element.node_i.x )
    │  ├─ y:  $(element.node_i.y )
    │  └─ z:  $(element.node_i.z )
    ├─ {bright_cyan:Node (j):}
    │  ├─ ID: $(element.node_j.ID)
    │  ├─ x:  $(element.node_j.x )
    │  ├─ y:  $(element.node_j.y )
    │  └─ z:  $(element.node_j.z )
    ├─ {bright_cyan:Section:}
    │  ├─ ID:   $(element.section.ID  )
    │  ├─ A:    $(element.section.A   )
    │  ├─ I_zz: $(element.section.I_zz)
    │  ├─ I_yy: $(element.section.I_yy)
    │  └─ J:    $(element.section.J   )
    ├─ {bright_cyan:Material:}
    │  ├─ ID: $(element.material.ID)
    │  ├─ E:  $(element.material.E )
    │  ├─ ν:  $(element.material.ν )
    │  └─ ρ:  $(element.material.ρ )
    └─ ω: $(element.ω)°""")

Base.show(io::IO, concload::ConcentratedLoad) = print(io, 
    styled"{cyan, bold:Concentrated load} applied to node with ID $(concload.ID)")
Base.show(io::IO, ::MIME"text/plain", concload::ConcentratedLoad) = print(io,
    styled"""
    {cyan, bold:Concentrated load:}
    ├─ Applied to node with ID: $(concload.ID)
    ├─ Forces:
    │  ├─ F_x: $(concload.F_x)
    │  ├─ F_y: $(concload.F_y)
    │  └─ F_z: $(concload.F_z)
    └─ Moments:
       ├─ M_x: $(concload.M_x)
       ├─ M_y: $(concload.M_y)
       └─ M_z: $(concload.M_z)""")

Base.show(io::IO, distload::DistributedLoad) = print(io,
    styled"{cyan, bold:Distributed load} applied to element with ID $(distload.ID)")
Base.show(io::IO, ::MIME"text/plain", distload::DistributedLoad) = print(io,
    styled"""
    {cyan, bold:Distributed load:}
    ├─ Applied to element with ID: $(distload.ID)
    ├─ w_x: $(distload.w_x)
    ├─ w_y: $(distload.w_y)
    └─ w_z: $(distload.w_z)""")

function Base.show(io::IO, ::MIME"text/plain", model::Model)
    emptyflag =     
        isempty(model.nodes    ) && 
        isempty(model.sections ) && 
        isempty(model.materials) && 
        isempty(model.elements ) && 
        isempty(model.concloads) && 
        isempty(model.distloads)

    if emptyflag
        print(io, styled"{yellow, bold:Empty model.}")
    else
        numnodes     = length(model.nodes    )
        numsections  = length(model.sections )
        nummaterials = length(model.materials)
        numelements  = length(model.elements )
        numconcloads = length(model.concloads)
        numdistloads = length(model.distloads)

        print(io, styled"""
        {cyan, bold:Model:}
        ├─ Nodes:     $(numnodes    )
        ├─ Sections:  $(numsections )
        ├─ Materials: $(nummaterials)
        ├─ Elements:  $(numelements )
        └─ Loads:
           ├─ Concentrated: $(numconcloads)
           └─ Distributed:  $(numdistloads)""")
    end
end