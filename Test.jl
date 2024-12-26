# Reference:
# http://www1.coe.neu.edu/~jfhajjar/home/Denavit%20and%20Hajjar%20-%20Geometric%20Nonlinearity%20in%20OpenSees%20-%20Report%20No.%20NEU-CEE-2013-02%202013.pdf

using Hephaestus
using ReverseDiff

function f(x)
# Define an empty model:
M = Model()

# Define the nodes and DOF supports:
node!(M,  1,   0, 0, 0, u_x = true, u_y = true, u_z = true, θ_x = true, θ_y = true)
node!(M,  2,  12, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  3,  24, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  4,  36, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  5,  48, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  6,  60, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  7,  72, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  8,  84, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M,  9,  96, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M, 10, 108, 0, 0, u_z = true, θ_x = true, θ_y = true)
node!(M, 11, 120, 0, 0, u_y = true, u_z = true, θ_x = true, θ_y = true)

# Define the sections:
section!(M, 1, 10, 100, 100, 5)

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
distload!(M,  1, 0, -x[1], 0)
distload!(M,  2, 0, -x[1], 0)
distload!(M,  3, 0, -x[1], 0)
distload!(M,  4, 0, -x[1], 0)
distload!(M,  5, 0, -x[1], 0)
distload!(M,  6, 0, -x[1], 0)
distload!(M,  7, 0, -x[1], 0)
distload!(M,  8, 0, -x[1], 0)
distload!(M,  9, 0, -x[1], 0)
distload!(M, 10, 0, -x[1], 0)

U = solve(M, LinearElasticAnalysis()).U

return Δ = extract_node_disp(M, U, 6)[2]
end