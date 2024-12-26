# Reference:
# http://www1.coe.neu.edu/~jfhajjar/home/Denavit%20and%20Hajjar%20-%20Geometric%20Nonlinearity%20in%20OpenSees%20-%20Report%20No.%20NEU-CEE-2013-02%202013.pdf

using Hephaestus
using GLMakie

# Define an empty model:
M = Model()

# Define the nodes and DOF supports:
node!(M,  1,   0, 0, 0, u_x = true, u_y = true, u_z = true, θ_x = true, θ_y = true, θ_z = true)
node!(M,  2,  18, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  3,  36, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  4,  54, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  5,  72, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  6,  90, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  7, 108, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  8, 126, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  9, 144, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M, 10, 162, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M, 11, 180, 0, 0, u_z = true, θ_x = true, θ_y = true)

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

# Define the loads:
P = 50
H = 1
concload!(M, 11, -P, -H, 0, 0, 0, 0)

begin
    F = Figure()

    A = Axis3(F[1, 1])

    plotmodel!(A, M)

    display(F)
end