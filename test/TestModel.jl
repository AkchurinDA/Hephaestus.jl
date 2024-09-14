using Hephaestus

# Create a new model:
M = Model()

# Add nodes to the model:
add_node!(M, 1, 0, 1 * 12, 0)
add_node!(M, 2, 0, 2 * 12, 0)
add_node!(M, 3, 0, 3 * 12, 0)
add_node!(M, 4, 0, 4 * 12, 0)
add_node!(M, 5, 0, 5 * 12, 0)

# Add fixities to the model:
add_fixity!(M, 1, true, true, true, true, true, true)

# Add materials to the model:
add_material!(M, 1, 29000, 0.3)

# Add sections to the model:
add_section!(M, 1, 10, 100, 100, 15)

# Add elements to the model:
add_element!(M, 1, 1, 2, 1, 1)
add_element!(M, 2, 2, 3, 1, 1)
add_element!(M, 3, 3, 4, 1, 1)
add_element!(M, 4, 4, 5, 1, 1)