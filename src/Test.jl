# Reference:
# http://www1.coe.neu.edu/~jfhajjar/home/Denavit%20and%20Hajjar%20-%20Geometric%20Nonlinearity%20in%20OpenSees%20-%20Report%20No.%20NEU-CEE-2013-02%202013.pdf

using Hephaestus
using ForwardDiff

function f(x)
    # Define an empty model:
    M = Model()

    # Define the nodes and DOF supports:
    node!(M,  1,   0, 0, 0, u_x = true, u_y = true, u_z = true, θ_x = true, θ_y = true, θ_z = true)
    node!(M,  2,  12, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M,  3,  24, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M,  4,  36, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M,  5,  48, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M,  6,  60, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M,  7,  72, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M,  8,  84, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M,  9,  96, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M, 10, 108, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M, 11, 120, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M, 12, 132, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M, 13, 144, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M, 14, 156, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M, 15, 168, 0, 0, u_z = true, θ_x = true, θ_y = true)
    node!(M, 16, 180, 0, 0, u_z = true, θ_x = true, θ_y = true)

    # Define the sections:
    section!(M, 1, 9.12, 110, 37.1, 0)

    # Define the materials:
    material!(M, 1, 29000, 0.3, 0)

    # Define the elements:
    element!(M,  1,  1,  2, 1, 1)
    element!(M,  2,  2,  3, 1, 1)
    element!(M,  3,  3,  4, 1, 1)
    element!(M,  4,  4,  5, 1, 1)
    element!(M,  5,  5,  6, 1, 1)
    element!(M,  6,  6,  7, 1, 1)
    element!(M,  7,  7,  8, 1, 1)
    element!(M,  8,  8,  9, 1, 1)
    element!(M,  9,  9, 10, 1, 1)
    element!(M, 10, 10, 11, 1, 1)
    element!(M, 11, 11, 12, 1, 1)
    element!(M, 12, 12, 13, 1, 1)
    element!(M, 13, 13, 14, 1, 1)
    element!(M, 14, 14, 15, 1, 1)
    element!(M, 15, 15, 16, 1, 1)

    # Define the loads:
    P = x
    H = 1
    concload!(M, 16, -P, -H, 0, 0, 0, 0)

    # Solve the model using a linear elastic analysis:
    U = solve(M, LinearElasticAnalysis())

    # Extract the vertical displacement of free end of the cantilever beam:
    Δ = extract_node_disp(M, U, 16)[2]
end

ForwardDiff.derivative(f, 50)