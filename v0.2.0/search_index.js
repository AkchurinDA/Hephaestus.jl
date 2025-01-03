var documenterSearchIndex = {"docs":
[{"location":"SectionLibrary/#Section-library","page":"Section library","title":"Section library","text":"","category":"section"},{"location":"SectionLibrary/#Generic-section","page":"Section library","title":"Generic section","text":"","category":"section"},{"location":"SectionLibrary/#Rectangular-section","page":"Section library","title":"Rectangular section","text":"","category":"section"},{"location":"SectionLibrary/#Circular-section","page":"Section library","title":"Circular section","text":"","category":"section"},{"location":"SectionLibrary/#I-section","page":"Section library","title":"I-section","text":"","category":"section"},{"location":"SectionLibrary/#T-section","page":"Section library","title":"T-section","text":"","category":"section"},{"location":"SectionLibrary/#L-section","page":"Section library","title":"L-section","text":"","category":"section"},{"location":"ElementLibrary/#Element-Library","page":"Element library","title":"Element Library","text":"","category":"section"},{"location":"ElementLibrary/#Euler-Bernoulli-beam-column-element","page":"Element library","title":"Euler-Bernoulli beam-column element","text":"","category":"section"},{"location":"ElementLibrary/#Timoshenko-beam-column-element","page":"Element library","title":"Timoshenko beam-column element","text":"","category":"section"},{"location":"ElementLibrary/#Truss-element","page":"Element library","title":"Truss element","text":"","category":"section"},{"location":"Tutorial/#Postprocessing","page":"Postprocessing","title":"Postprocessing","text":"","category":"section"},{"location":"Tutorial/#Extracting-nodal-displacements","page":"Postprocessing","title":"Extracting nodal displacements","text":"","category":"section"},{"location":"Tutorial/#Extracting-nodal-reactions","page":"Postprocessing","title":"Extracting nodal reactions","text":"","category":"section"},{"location":"Tutorial/#Extracting-element-displacements","page":"Postprocessing","title":"Extracting element displacements","text":"","category":"section"},{"location":"Tutorial/#Extract-element-forces","page":"Postprocessing","title":"Extract element forces","text":"","category":"section"},{"location":"Tutorial/#Generating-reports","page":"Postprocessing","title":"Generating reports","text":"","category":"section"},{"location":"MaterialLibrary/#Material-library","page":"Material library","title":"Material library","text":"","category":"section"},{"location":"MaterialLibrary/#Elastic-material","page":"Material library","title":"Elastic material","text":"","category":"section"},{"location":"MaterialLibrary/#Elastic-perfectly-plastic-material","page":"Material library","title":"Elastic-perfectly plastic material","text":"","category":"section"},{"location":"QuickStart/#Quick-Start","page":"Quick Start","title":"Quick Start","text":"","category":"section"},{"location":"QuickStart/#Defining-a-model","page":"Quick Start","title":"Defining a model","text":"","category":"section"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To quickly get you started with Hephaestus.jl, let us recreate a simple example of a cantilever beam subjected to a distributed load as shown below.","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"First of all, load the package using the following command:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"using Hephaestus","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To create a new model, use the Model() constructor:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"model = Model()","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To add nodes to the model, use the node!() function:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"node!(model,  1,  0 * 12, 0, 0, u_x = true, u_y = true, θ_z = true)\nnode!(model,  2,  1 * 12, 0, 0)\nnode!(model,  3,  2 * 12, 0, 0)\nnode!(model,  4,  3 * 12, 0, 0)\nnode!(model,  5,  4 * 12, 0, 0)\nnode!(model,  6,  5 * 12, 0, 0)\nnode!(model,  7,  6 * 12, 0, 0)\nnode!(model,  8,  7 * 12, 0, 0)\nnode!(model,  9,  8 * 12, 0, 0)\nnode!(model, 10,  9 * 12, 0, 0)\nnode!(model, 11, 10 * 12, 0, 0)","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To add sections to the model, use the section!() function:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"section!(model, 1, 9.16, 180, 180, 359)","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To add materials to the model, use the material!() function:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"material!(model, 1, 29000, 0.3, 0.290)","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To add elements to the model, use the element!() function:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"element!(model, 1 , 1 , 2 , 1, 1)\nelement!(model, 2 , 2 , 3 , 1, 1)\nelement!(model, 3 , 3 , 4 , 1, 1)\nelement!(model, 4 , 4 , 5 , 1, 1)\nelement!(model, 5 , 5 , 6 , 1, 1)\nelement!(model, 6 , 6 , 7 , 1, 1)\nelement!(model, 7 , 7 , 8 , 1, 1)\nelement!(model, 8 , 8 , 9 , 1, 1)\nelement!(model, 9 , 9 , 10, 1, 1)\nelement!(model, 10, 10, 11, 1, 1)","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"To add distributed loads to the model, use the distload!() function:","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"distload!(model,  1, 0, -1, 0)\ndistload!(model,  2, 0, -1, 0)\ndistload!(model,  3, 0, -1, 0)\ndistload!(model,  4, 0, -1, 0)\ndistload!(model,  5, 0, -1, 0)\ndistload!(model,  6, 0, -1, 0)\ndistload!(model,  7, 0, -1, 0)\ndistload!(model,  8, 0, -1, 0)\ndistload!(model,  9, 0, -1, 0)\ndistload!(model, 10, 0, -1, 0)","category":"page"},{"location":"QuickStart/#Plotting-model-and-solution","page":"Quick Start","title":"Plotting model and solution","text":"","category":"section"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"using CairoMakie","category":"page"},{"location":"QuickStart/","page":"Quick Start","title":"Quick Start","text":"plotmodel(model)","category":"page"},{"location":"QuickStart/#Generating-reports","page":"Quick Start","title":"Generating reports","text":"","category":"section"},{"location":"#Description","page":"Home","title":"Description","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Hephaestus.jl is an automatically differentiable structural analysis package purely written in the Julia programming language.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install Hephaestus.jl package, type ] in Julia REPL to enter the built-in Julia package manager and execute the following command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add Hephaestus","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Hephaestus.jl package is distributed under the MIT license. More information can be found in the LICENSE.md file.","category":"page"},{"location":"#Help-and-Support","page":"Home","title":"Help and Support","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For assistance with the package, please raise an issue on the GitHub Issues page. Please use the appropriate labels to indicate the specific functionality you are inquiring about. Alternatively, contact the author directly at AkchurinDA@gmail.com.","category":"page"},{"location":"#Acknowledgements","page":"Home","title":"Acknowledgements","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The internal design of the package is inspired by OpenSeesPy, PyNite, and MASTAN2.","category":"page"},{"location":"API/#API","page":"API","title":"API","text":"","category":"section"},{"location":"API/#Types-Used-to-Define-a-Model","page":"API","title":"Types Used to Define a Model","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Model\nNode\nMaterial\nSection\nElement","category":"page"},{"location":"API/#Hephaestus.Model","page":"API","title":"Hephaestus.Model","text":"struct Model\n\nA type representing the finite element model of a structure of interest.\n\ndimensionality: Dimensionality of the model\nnodes: Nodes of the model\nsections: Sections of the model\nmaterials: Materials of the model\nelements: Elements of the model\nconcloads: Concentrated loads of the model\ndistloads: Distribution loads of the model\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.Node","page":"API","title":"Hephaestus.Node","text":"struct Node\n\nA type representing a node in a finite element model.\n\nID: Identification tag\nx: x-coordinate\ny: y-coordinate\nz: z-coordinate\nu_x: Is the node restrained against translation along the global x-axis?\nu_y: Is the node restrained against translation along the global y-axis?\nu_z: Is the node restrained against translation along the global z-axis?\nθ_x: Is the node restrained against rotation about the global x-axis?\nθ_y: Is the node restrained against rotation about the global y-axis?\nθ_z: Is the node restrained against rotation about the global z-axis?\nstate: Current state\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.Material","page":"API","title":"Hephaestus.Material","text":"struct Material\n\nA type representing a material in a finite element model.\n\nID: Identification tag\nE: Young's modulus, E\nν: Poisson's ratio, nu\nρ: Density, rho\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.Section","page":"API","title":"Hephaestus.Section","text":"struct Section\n\nA type representing a section in a finite element model.\n\nID: Identification tag\nA: Cross-sectional area, A\nI_zz: Moment of inertia about the local z-axis, I_zz\nI_yy: Moment of inertia about the local y-axis, I_yy\nJ: Torsional constant, J\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.Element","page":"API","title":"Hephaestus.Element","text":"struct Element\n\nA type representing an element in the finite element model.\n\nID: Identification tag\nnode_i: Node (i) of the element\nnode_j: Node (j) of the element\nsection: Section attached to the element\nmaterial: Material attached to the element\nω: Orientation angle of the element\nreleases_i: End moment releases at node (i) of the element\nreleases_j: End moment releases at node (j) of the element\nstate: Current state\n\n\n\n\n\n","category":"type"},{"location":"API/#Types-Used-to-Perform-Analyses-of-Different-Types-and-Store-the-Results","page":"API","title":"Types Used to Perform Analyses of Different Types and Store the Results","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"LinearElasticAnalysis\nLinearElasticAnalysisCache\nElasticBucklingAnalysis\nElasticBucklingAnalysisCache\nFreeVibrationAnalysis\nFreeVibrationAnalysisCache","category":"page"},{"location":"API/#Hephaestus.LinearElasticAnalysis","page":"API","title":"Hephaestus.LinearElasticAnalysis","text":"struct LinearElasticAnalysis\n\nA type representing the (geometrically) linear (materially) elastic analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.LinearElasticAnalysisCache","page":"API","title":"Hephaestus.LinearElasticAnalysisCache","text":"struct LinearElasticAnalysisCache\n\nA type used to store the results of linear elastic analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.ElasticBucklingAnalysis","page":"API","title":"Hephaestus.ElasticBucklingAnalysis","text":"struct ElasticBucklingAnalysis\n\nA type representing the elastic buckling analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.ElasticBucklingAnalysisCache","page":"API","title":"Hephaestus.ElasticBucklingAnalysisCache","text":"struct ElasticBucklingAnalysisCache\n\nA type used to store the results of elastic buckling analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.FreeVibrationAnalysis","page":"API","title":"Hephaestus.FreeVibrationAnalysis","text":"struct LinearElasticAnalysis\n\nA type representing the free vibration analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Hephaestus.FreeVibrationAnalysisCache","page":"API","title":"Hephaestus.FreeVibrationAnalysisCache","text":"struct LinearElasticAnalysisCache\n\nA type used to store the results of free vibration analysis.\n\n\n\n\n\n","category":"type"},{"location":"API/#Functions-Used-to-Define-a-Model","page":"API","title":"Functions Used to Define a Model","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"node!\nsection!\nmaterial!\nelement!\nconcload!\ndistload!","category":"page"},{"location":"API/#Hephaestus.node!","page":"API","title":"Hephaestus.node!","text":"node!(model::Model, ID::Int,\n    x::Real, y::Real, z::Real;\n    u_x::Bool = false, u_y::Bool = false, u_z::Bool = false,\n    θ_x::Bool = false, θ_y::Bool = false, θ_z::Bool = false)::Model\n\nAdd a node to a finite element model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.section!","page":"API","title":"Hephaestus.section!","text":"section!(model::Model, ID::Int,\n    A::Real, I_zz::Real, I_yy::Real, J::Real)::Model\n\nAdd a section to a finite element model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.material!","page":"API","title":"Hephaestus.material!","text":"material!(model::Model, ID::Int,\n    E::Real, ν::Real, ρ::Real)::Model\n\nAdd a material to a finite element model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.element!","page":"API","title":"Hephaestus.element!","text":"element!(model::Model, ID::Int,\n    node_i_ID  ::Int,\n    node_j_ID  ::Int,\n    section_ID ::Int,\n    material_ID::Int;\n    ω          ::Real = 0,\n    releases_i ::Vector{<:Bool} = [false, false, false],\n    releases_j ::Vector{<:Bool} = [false, false, false])::Model\n\nAdd an element to a finite element model.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.concload!","page":"API","title":"Hephaestus.concload!","text":"concload!(model::Model, ID::Int,\n    F_x::Real, F_y::Real, F_z::Real,\n    M_x::Real, M_y::Real, M_z::Real)::Model\n\nApplies a concentrated load to a node with a specified ID.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.distload!","page":"API","title":"Hephaestus.distload!","text":"distload!(model::Model, ID::Int,\n    w_x::Real, w_y::Real, w_z::Real)::Model\n\nApplies a distributed load to an element with a specified ID.\n\n\n\n\n\n","category":"function"},{"location":"API/#Functions-Used-to-Perform-Analyses-of-Different-Types-and-Extract-the-Results","page":"API","title":"Functions Used to Perform Analyses of Different Types and Extract the Results","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"solve!\ngetnodaldisplacements\ngetnodalreactions\ngetelementdisplacements\ngetelementforces","category":"page"},{"location":"API/#Hephaestus.solve!","page":"API","title":"Hephaestus.solve!","text":"solve!(model::Model, analysistype::AbstractAnalysisType; continueanalaysis::Bool = false)\n\nSolve the model using the specified analysis type.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.getnodaldisplacements","page":"API","title":"Hephaestus.getnodaldisplacements","text":"getnodaldisplacements(model::Model, ID::Int)\n\nExtracts the displacement vector of a node of interest.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.getnodalreactions","page":"API","title":"Hephaestus.getnodalreactions","text":"getnodalreactions(model::Model, ID::Int)\n\nExtracts the reaction vector of a node of interest from the solution cache.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.getelementdisplacements","page":"API","title":"Hephaestus.getelementdisplacements","text":"getelementdisplacements(model::Model, ID::Int)\n\nExtracts the displacement vector of an element of interest in its local coordinate system.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.getelementforces","page":"API","title":"Hephaestus.getelementforces","text":"getelementforces(model::Model, ID::Int)\n\nExtracts the internal force vector of an element of interest in its local coordinate system.\n\n\n\n\n\n","category":"function"},{"location":"API/#Functions-Used-to-Plot-a-Model-and-the-Results-of-Analyses-of-Different-Types","page":"API","title":"Functions Used to Plot a Model and the Results of Analyses of Different Types","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"plotmodel\nplotmodel!","category":"page"},{"location":"API/#Hephaestus.plotmodel","page":"API","title":"Hephaestus.plotmodel","text":"plotmodel(model::Model, [options])\n\nPlots a model of a structure of interest (Model) into a new Makie.jl scene.\n\n\n\n\n\n","category":"function"},{"location":"API/#Hephaestus.plotmodel!","page":"API","title":"Hephaestus.plotmodel!","text":"plotmodel!(scene::Scene, model::Model, [options])\n\nPlots a model of a structure of interest (Model) into an existing Makie.jl scene.\n\n\n\n\n\n","category":"function"}]
}