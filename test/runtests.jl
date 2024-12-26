using Test
using Hephaestus
using DifferentiationInterface, ForwardDiff, ReverseDiff

@testset "Linear Elastic Analysis" begin
    # Reference:
    # http://www1.coe.neu.edu/~jfhajjar/home/Denavit%20and%20Hajjar%20-%20Geometric%20Nonlinearity%20in%20OpenSees%20-%20Report%20No.%20NEU-CEE-2013-02%202013.pdf
    
    # Define an empty model:
    model = Model()

    # Define the nodes and DOF supports:
    node!(model,  1,   0, 0, 0, u_x = true, u_y = true, u_z = true, θ_x = true, θ_y = true, θ_z = true)
    node!(model,  2,  12, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model,  3,  24, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model,  4,  36, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model,  5,  48, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model,  6,  60, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model,  7,  72, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model,  8,  84, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model,  9,  96, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model, 10, 108, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model, 11, 120, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model, 12, 132, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model, 13, 144, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model, 14, 156, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model, 15, 168, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(model, 16, 180, 0, 0, u_z = true, θ_x = true, θ_y = true)

    # Define the sections:
    section!(model, 1, 9.12, 110, 37.1, 0)

    # Define the materials:
    material!(model, 1, 29000, 0.3, 0)

    # Define the elements:
    element!(model,  1,  1,  2, 1, 1)
    element!(model,  2,  2,  3, 1, 1)
    element!(model,  3,  3,  4, 1, 1)
    element!(model,  4,  4,  5, 1, 1)
    element!(model,  5,  5,  6, 1, 1)
    element!(model,  6,  6,  7, 1, 1)
    element!(model,  7,  7,  8, 1, 1)
    element!(model,  8,  8,  9, 1, 1)
    element!(model,  9,  9, 10, 1, 1)
    element!(model, 10, 10, 11, 1, 1)
    element!(model, 11, 11, 12, 1, 1)
    element!(model, 12, 12, 13, 1, 1)
    element!(model, 13, 13, 14, 1, 1)
    element!(model, 14, 14, 15, 1, 1)
    element!(model, 15, 15, 16, 1, 1)

    # Define the loads:
    concload!(model, 16, -50, -1, 0, 0, 0, 0)

    # Solve the model using a linear elastic analysis:
    solution = solve(model, LinearElasticAnalysis())

    # Extract the vertical displacement of free end of the cantilever beam:
    Δ = getnodedisp(model, solution, 16)[2]

    # Check the result:
    @test Δ ≈ -0.609 atol = 1E-3
end
