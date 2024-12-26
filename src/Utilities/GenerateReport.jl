# What is included in the report:
# Header:
# ─ Date and time of the report
# ─ Name of the model
# Body:
# ─ Number of model entities

"""
    generatereport(model::Model; path = pwd(), reportname = "Report")

Generate a report.
"""
function generatereport(model::Model; path = pwd(), reportname = "Report")
    open(joinpath(path, "$(reportname).txt"), "w") do io
        write(io, "Hephaestus.jl: \n")
        write(io, "Automatically Differentiable \n")
        write(io, "Structural Analysis in Julia \n")
        write(io, "\n")

        currentdatetime = Dates.now()
        currentdate     = Dates.format(currentdatetime, "u. dd, yyyy")
        currenttime     = Dates.format(currentdatetime, "II:MM:SS p")
        write(io, "HEADER ******************************************* \n")
        write(io, "Model name: $(model.name) \n")
        write(io, "Date: $(currentdate) \n")
        write(io, "Time: $(currenttime) \n")
        write(io, "************************************************** \n")
        write(io, "\n")

        write(io, "MODEL ******************************************** \n")
        write(io, "GENERAL INFORMATION:                               \n")
        write(io, "(1) # of nodes:       $(length(model.nodes    ))   \n")
        write(io, "(2) # of sections:    $(length(model.sections ))   \n")
        write(io, "(3) # of materials:   $(length(model.materials))   \n")
        write(io, "(4) # of elements:    $(length(model.elements ))   \n")
        write(io, "(5) # of conc. loads: $(length(model.concloads))   \n")
        write(io, "(6) # of dist. loads: $(length(model.distloads))   \n")
        write(io, "\n")

        write(io, "(1) NODES: \n")
        write(io, "\n")
        if !isempty(model.nodes)
            write(io, "Nodal coordinates: \n")
            write(io, "│ ID         │ x          │ y          │ z          │ \n")
            write(io, "├────────────┼────────────┼────────────┼────────────┤ \n")
            for node in model.nodes
                @printf(io, "│ %10d │ %10.3e │ %10.3e │ %10.3e │ \n", 
                    node.ID, node.x, node.y, node.z)
            end
            write(io, "\n")
            write(io, "DOF supports: \n")
            write(io, "│ ID         │ u_x   │ u_y   │ u_z   │ θ_x   │ θ_y   │ θ_z   │ \n")
            write(io, "├────────────┼───────┼───────┼───────┼───────┼───────┼───────┤ \n")
            for node in model.nodes
                u_x_support = node.u_x == false ? "Free" : "Fixed"
                u_y_support = node.u_y == false ? "Free" : "Fixed"
                u_z_support = node.u_z == false ? "Free" : "Fixed"
                θ_x_support = node.θ_x == false ? "Free" : "Fixed"
                θ_y_support = node.θ_y == false ? "Free" : "Fixed"
                θ_z_support = node.θ_z == false ? "Free" : "Fixed"
                @printf(io, "│ %10d │ %5s │ %5s │ %5s │ %5s │ %5s │ %5s │ \n", 
                    node.ID, u_x_support, u_y_support, u_z_support, θ_x_support, θ_y_support, θ_z_support)
            end
        else
            write(io, "No nodes defined. \n")
        end
        write(io, "\n")
        write(io, "(2) SECTIONS: \n")
        write(io, "\n")
        if !isempty(model.sections)
            write(io, "│ ID         │ A          │ I_zz       │ I_yy       │ J          │ \n")
            write(io, "├────────────┼────────────┼────────────┼────────────┼────────────┤ \n")
            for section in model.sections
                @printf(io, "│ %10d │ %10.3e │ %10.3e │ %10.3e │ %10.3e │ \n", 
                    section.ID, section.A, section.I_zz, section.I_yy, section.J)
            end
        else
            write(io, "No sections defined. \n")
        end
        write(io, "\n")
        write(io, "(3) MATERIALS: \n")
        write(io, "\n")
        if !isempty(model.materials)
            write(io, "│ ID         │ E          │ ν          │ ρ          │ \n")
            write(io, "├────────────┼────────────┼────────────┼────────────┤ \n")
            for material in model.materials
                @printf(io, "│ %10d │ %10.3e │ %10.3e │ %10.3e │ \n", 
                    material.ID, material.E, material.ν, material.ρ)
            end
        else
            write(io, "No materials defined. \n")
        end
        write(io, "\n")
        write(io, "(4) ELEMENTS: \n")
        write(io, "\n")
        if !isempty(model.elements)
            write(io, "│ ID         │ Node (i)   │ Node (j)   │ Section    │ Material   │ ω (°)      │ \n")
            write(io, "├────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤ \n")
            for element in model.elements
                @printf(io, "│ %10d │ %10d │ %10d │ %10d │ %10d │ %10.3e │ \n", 
                    element.ID, element.node_i.ID, element.node_j.ID, element.section.ID, element.material.ID, element.ω)
            end
        else
            write(io, "No elements defined. \n")
        end
        write(io, "\n")
        write(io, "(5) CONC. LOADS: \n")
        write(io, "\n")
        if !isempty(model.concloads)
            write(io, "│ Node ID    │ F_x        │ F_y        │ F_z        │ M_x        │ M_y        │ M_z        │ \n")
            write(io, "├────────────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤ \n")
            for concload in model.concloads
                @printf(io, "│ %10d │ %10.3e │ %10.3e │ %10.3e │ %10.3e │ %10.3e │ %10.3e │ \n", 
                    concload.ID, concload.F_x, concload.F_y, concload.F_z, concload.M_x, concload.M_y, concload.M_z)
            end
        else
            write(io, "No concentrated loads defined. \n")
        end
        write(io, "\n")
        write(io, "(6) DIST. LOADS: \n")
        write(io, "\n")
        if !isempty(distload!)
            write(io, "│ Element ID │ w_x        │ w_y        │ w_z        │ CS      │ \n")
            write(io, "├────────────┼────────────┼────────────┼────────────┼─────────┤ \n")
            for distload in model.distloads
                cs = distload.cs == :global ? "Global" : "Local"
                @printf(io, "│ %10d │ %10.3e │ %10.3e │ %10.3e │ %7s │ \n",
                    distload.ID, distload.w_x, distload.w_y, distload.w_z, cs)
            end
        else
            write(io, "No distributed loads defined. \n")
        end
        write(io, "************************************************** \n")
        write(io, "\n")
    end

    return nothing
end
