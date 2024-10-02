var documenterSearchIndex = {"docs":
[{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To quickly get you started with Hephaestus.jl, let us recreate a simple example of a cantilever beam subjected to a concentrated load at its free end as shown below.","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"using Hephaestus","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To create a new model, use the Model() constructor:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"M = Model()","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To add nodes to the model, use the add_node!() function:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"add_node!(M, 1 , 0  * 12, 0, 0)\nadd_node!(M, 2 , 1  * 12, 0, 0)\nadd_node!(M, 3 , 2  * 12, 0, 0)\nadd_node!(M, 4 , 3  * 12, 0, 0)\nadd_node!(M, 5 , 4  * 12, 0, 0)\nadd_node!(M, 6 , 5  * 12, 0, 0)\nadd_node!(M, 7 , 6  * 12, 0, 0)\nadd_node!(M, 8 , 7  * 12, 0, 0)\nadd_node!(M, 9 , 8  * 12, 0, 0)\nadd_node!(M, 10, 9  * 12, 0, 0)\nadd_node!(M, 11, 10 * 12, 0, 0)","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To add materials to the model, use the add_material!() function:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"add_material!(M, 1, 29000, 0.3, 0.000290)","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To add sections to the model, use the add_section!() function:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"add_section!(M, 1, 9.16, 180, 180, 359)","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To add elements to the model, use the add_element!() function:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"add_element!(M, 1 , 1 , 2 , 1, 1)\nadd_element!(M, 2 , 2 , 3 , 1, 1)\nadd_element!(M, 3 , 3 , 4 , 1, 1)\nadd_element!(M, 4 , 4 , 5 , 1, 1)\nadd_element!(M, 5 , 5 , 6 , 1, 1)\nadd_element!(M, 6 , 6 , 7 , 1, 1)\nadd_element!(M, 7 , 7 , 8 , 1, 1)\nadd_element!(M, 8 , 8 , 9 , 1, 1)\nadd_element!(M, 9 , 9 , 10, 1, 1)\nadd_element!(M, 10, 10, 11, 1, 1)","category":"page"},{"location":"","page":"Home","title":"Home","text":"<img src = \"./assets/social-preview.svg\"/>","category":"page"},{"location":"#Description","page":"Home","title":"Description","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Hephaestus.jl is an auto-differentiable structural analysis package purely written in the Julia programming language.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install Hephaestus.jl package, type ] in Julia REPL to enter the built-in Julia package manager and execute the following command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add Hephaestus","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Hephaestus.jl package is distributed under the MIT license. More information can be found in the LICENSE.md file.","category":"page"},{"location":"#Help-and-Support","page":"Home","title":"Help and Support","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For assistance with the package, please raise an issue on the GitHub Issues page. Please use the appropriate labels to indicate the specific functionality you are inquiring about. Alternatively, contact the author directly at AkchurinDA@gmail.com.","category":"page"},{"location":"#Acknowledgements","page":"Home","title":"Acknowledgements","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The design of the package is inspired by OpenSeesPy, PyNite, and MASTAN2.","category":"page"},{"location":"API/#Types","page":"API","title":"Types","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Model\nNode\nMaterial\nSection\nElement\nO1EAnalysis\nO1ESolutionCache\nO2EAnalysis\nO2ESolutionCache\nEBAnalysis\nEBSolutionCache","category":"page"},{"location":"API/#Hephaestus.Model","page":"API","title":"Hephaestus.Model","text":"struct Model\n\nA type representing the finite element model of a structure of interest.\n\nnodes: Ordered dictionary that stores the nodes of the model\nmaterials: Ordered dictionary that stores the materials of the model\nsections: Ordered dictionary that stores the sections of the model\nelements: Ordered dictionary that stores the elements of the model\nsprings: Ordered dictionary that stores the springs of the model\nsupports: Ordered dictionary that stores the supports of the model\nconc_loads: Ordered dictionary that stores the concentrated loads of the model\ndist_loads: Ordered dictionary that stores the distributed loads of the model\np_l: Ordered dictionary that stores the element fixed-end forces in the local coordinate system\np_g: Ordered dictionary that stores the element fixed-end forces in the global coordinate system\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.Node","page":"API","title":"Hephaestus.Node","text":"struct Node\n\nA type representing a node in the model of a structure of interest.\n\nID: Unique identifier\nx: x-coordinate\ny: y-coordinate\nz: z-coordinate\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.Material","page":"API","title":"Hephaestus.Material","text":"struct Material\n\nA type representing a material in the model of a structure of interest.\n\nID: Unique identifier\nE: Young's modulus, E\nν: Poisson's ratio, u\nρ: Density, rho\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.Section","page":"API","title":"Hephaestus.Section","text":"struct Section\n\nA type representing a section in the model of a structure of interest.\n\nID: Unique identifier\nA: Cross-sectional area, A\nI_zz: Moment of inertia about the local z-axis, I_zz\nI_yy: Moment of inertia about the local y-axis, I_yy\nJ: Polar moment of inertia about the local x-axis, J\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.Element","page":"API","title":"Hephaestus.Element","text":"struct Element\n\nA type representing an element in the model of a structure of interest.\n\nID: Unique identifier\nnode_i_ID: Unique identifier of the node i of the element\nnode_j_ID: Unique identifier of the node j of the element\nmaterial_ID: Unique identifier of the material of the element\nsection_ID: Unique identifier of the section of the element\nreleases_i: DOF releases at the node i\nreleases_j: DOF releases at the node j\nω: Angle that defines the orientation of the element's local coordinate system, omega\nx_i: x-coordinate of the node i, x_i\ny_i: y-coordinate of the node i, y_i\nz_i: z-coordinate of the node i, z_i\nx_j: x-coordinate of the node j, x_j\ny_j: y-coordinate of the node j, y_j\nz_j: z-coordinate of the node j, z_j\nE: Young's modulus, E\nν: Poisson's ratio, nu\nρ: Density, rho\nA: Cross-sectional area, A\nI_zz: Moment of inertia about the local z-axis, I_zz\nI_yy: Moment of inertia about the local y-axis, I_yy\nJ: Polar moment of inertia about the local x-axis, J\nL: Length of the element, L\nγ: Local-to-global sub-transformation matrix, gamma\nΓ: Local-to-global transformation matrix, Gamma\nk_e_l: Elastic stiffness matrix in the local coordinate system of the element, k_e l\nk_e_g: Elastic stiffness matrix in the global coordinate system, k_e g\nk_g_l: Geometric stiffness matrix in the local coordinate system of the element, k_g l\nk_g_g: Geometric stiffness matrix in the global coordinate system, k_g g\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.O1EAnalysis","page":"API","title":"Hephaestus.O1EAnalysis","text":"struct O1EAnalysis\n\nA type that represents the 1st-order (O1) elastic (E) analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.O1ESolutionCache","page":"API","title":"Hephaestus.O1ESolutionCache","text":"struct O1ESolutionCache\n\nA type that stores the results of the 1st-order elastic analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.O2EAnalysis","page":"API","title":"Hephaestus.O2EAnalysis","text":"struct O2EAnalysis\n\nA type that represents the 2nd-order (O2) elastic (E) analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.O2ESolutionCache","page":"API","title":"Hephaestus.O2ESolutionCache","text":"struct O2ESolutionCache\n\nA type that stores the results of the 2nd-order elastic analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.EBAnalysis","page":"API","title":"Hephaestus.EBAnalysis","text":"struct EBAnalysis\n\nA type that represents the elastic buckling (EB) analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.EBSolutionCache","page":"API","title":"Hephaestus.EBSolutionCache","text":"struct EBSolutionCache\n\nA type that stores the results of the elastic buckling analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Functions","page":"API","title":"Functions","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"add_node!\ndel_node!\nadd_material!\ndel_material!\nadd_section!\ndel_section!\nadd_element!\ndel_element!\nadd_support!\ndel_support!\nadd_conc_load!\ndel_conc_load!\nadd_dist_load!\ndel_dist_load!\nget_node_u_g\nget_element_u_l\nget_element_f_l\nplotmodel\nplotmodel!","category":"page"},{"location":"API/#Hephaestus.add_node!","page":"API","title":"Hephaestus.add_node!","text":"add_node!(model::Model, ID::Int, \n    x::Real, y::Real, z::Real)\n\nAdds a node to the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.del_node!","page":"API","title":"Hephaestus.del_node!","text":"del_node!(model::Model, ID::Int)\n\nDeletes a node from the model and all related elements, supports, and concentrated loads.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.add_material!","page":"API","title":"Hephaestus.add_material!","text":"add_material!(model::Model, ID::Int, \n    E::Real, ν::Real, ρ::Real)\n\nAdds a material to the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.del_material!","page":"API","title":"Hephaestus.del_material!","text":"del_material!(model::Model, ID::Int)\n\nDeletes a material from the model and all related elements.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.add_section!","page":"API","title":"Hephaestus.add_section!","text":"add_section!(model::Model, ID::Int, \n    A::Real, I_zz::Real, I_yy::Real, J::Real)\n\nAdds a section to the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.del_section!","page":"API","title":"Hephaestus.del_section!","text":"del_section!(model::Model, ID::Int)\n\nDeletes a section from the model and all related elements.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.add_element!","page":"API","title":"Hephaestus.add_element!","text":"add_element!(model::Model, ID::Int,\n    node_i_ID::Int, node_j_ID::Int, material_ID::Int, section_ID::Int;\n    releases_i::Vector{Bool} = [false, false, false, false, false, false],\n    releases_j::Vector{Bool} = [false, false, false, false, false, false],\n    ω::Real = 0)\n\nAdds an element to the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.del_element!","page":"API","title":"Hephaestus.del_element!","text":"del_element!(model::Model, ID::Int)\n\nDeletes an element from the model and all related distributed loads.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.add_support!","page":"API","title":"Hephaestus.add_support!","text":"add_support!(model::Model, ID::Int, \n    u_x::Bool, u_y::Bool, u_z::Bool, \n    θ_x::Bool, θ_y::Bool, θ_z::Bool)\n\nAdds a support to the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.del_support!","page":"API","title":"Hephaestus.del_support!","text":"del_support!(model::Model, ID::Int)\n\nDeletes a support from the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.add_conc_load!","page":"API","title":"Hephaestus.add_conc_load!","text":"add_conc_load!(model::Model, ID::Int,\n    F_x::Real, F_y::Real, F_z::Real,\n    M_x::Real, M_y::Real, M_z::Real)\n\nAdds a concentrated load to the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.del_conc_load!","page":"API","title":"Hephaestus.del_conc_load!","text":"del_conc_load!(model::Model, ID::Int)\n\nDeletes a concentrated load from the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.add_dist_load!","page":"API","title":"Hephaestus.add_dist_load!","text":"add_dist_load!(model::Model, ID::Int,\n    q_x::Real, q_y::Real, q_z::Real)\n\nAdds a distributed load to the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.del_dist_load!","page":"API","title":"Hephaestus.del_dist_load!","text":"del_dist_load!(model::Model, ID::Int)\n\nDeletes a distributed load from the model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.get_node_u_g","page":"API","title":"Hephaestus.get_node_u_g","text":"get_node_u_g(solution_cache::AbstractSolutionCache, ID::Int)\n\nExtracts the displacement vector of a node in the global coordinate system.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.get_element_u_l","page":"API","title":"Hephaestus.get_element_u_l","text":"get_element_u_l(model::Model, solution_cache::AbstractSolutionCache, ID::Int)\n\nExtracts the element displacement vector in the local coordinate system.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.get_element_f_l","page":"API","title":"Hephaestus.get_element_f_l","text":"get_element_f_l(model::Model, solution_cache::AbstractSolutionCache, ID::Int)\n\nExtract the element force vector in the local coordinate system.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.plotmodel","page":"API","title":"Hephaestus.plotmodel","text":"plotmodel(model::Model, [options])\n\nPlots a model of a structure of interest (Model) into a new Makie.jl scene.\n\nPlotting Options\n\nNodes\n\nOption Description Default\nnode_color Color of the nodes :red\nnode_marker Marker of the nodes :circle\nnode_markersize Size of the nodes 6\nnode_strokecolor Stroke color of the nodes :black\nnode_strokewidth Stroke width of the nodes 1\n\nNode Labels\n\nOption Description Default\nlabel_nodes Whether to label the nodes true\nalign_node_labels Alignment of the node labels (:center, :bottom)\nnode_label_color Color of the node labels :red\n\nElements\n\nOption Description Default\nelement_color Color of the elements :black\nelement_linestyle Line style of the elements :solid\nelement_linewidth Line width of the elements 1\n\nElement Labels\n\nOption Description Default\nlabel_elements Whether to label the elements true\nalign_element_labels Alignment of the element labels (:center, :bottom)\nelement_label_color Color of the element labels :black\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.plotmodel!","page":"API","title":"Hephaestus.plotmodel!","text":"plotmodel!(model::Model, [options])\n\nPlots a model of a structure of interest (Model) into an existing scene. The plotting options are the same as in plotmodel() function.\n\n\n\n\n\n","category":"function"}]
}
