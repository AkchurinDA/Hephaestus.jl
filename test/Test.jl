using Hephaestus

# Create a new model:
M = Model()

# Add nodes:
add_node!(M, 1, 0    , 0   , 0)
add_node!(M, 2, 0    , 5000, 0)
add_node!(M, 3, 7416 , 8000, 0)
add_node!(M, 4, 11416, 5000, 0)
add_node!(M, 5, 11416, 0   , 0)

# Add materials:
add_material!(M, 1, 200000, 0.3, 9.81 * 7850)

# Add sections:
add_section!(M, 1, 4E3, 50E6 , 0, 0)
add_section!(M, 2, 6E3, 200E6, 0, 0)

# Add elements:
add_element!(M, 1, 1, 2, 1, 1)
add_element!(M, 2, 2, 3, 1, 2)
add_element!(M, 3, 3, 4, 1, 1)
add_element!(M, 4, 4, 5, 1, 1)

add_nodal_load!(M, 2, 0, -5000.0, π, 0, 0, 0)

# Add supports:
add_support!(M, 1, true, true, true, true, true, true)

solve(M, :O1E)