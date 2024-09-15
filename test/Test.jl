using Hephaestus

M = Model()

add_node!(M, 1, 0.0, 0.0, 0.0)
add_node!(M, 2, 1.0, 0.0, 0.0)

add_material!(M, 1, 200000, 0.3, 100)

add_section!(M, 1, 8000, 200000000, 0, 300000)

add_element!(M, 1, 1, 2, 1, 1, ω = π / 4)
