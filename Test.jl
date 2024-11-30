using Hephaestus

model = Model()

add_node!(model, 1, 0, 0, 0)
add_node!(model, 2, 1, 0, 0)
add_node!(model, 3, 2, 0, 0)
add_node!(model, 4, 3, 0, 0)
add_node!(model, 5, 4, 0, 0)
add_node!(model, 100, 5, 0, 0)

add_material!(model, 1, 29000, 0.3, 0.290)

add_section!(model, 1, 10, 100, 100, 5)

add_element!(model, 1, 1, 2, 1, 1)
add_element!(model, 2, 2, 3, 1, 1)
add_element!(model, 3, 3, 4, 1, 1)
add_element!(model, 4, 4, 5, 1, 1)
add_element!(model, 5, 5, 100, 1, 1)

add_support!(model, 1, true, true, true, true, true, true)

add_cload!(model, 100, 0, -10, 0, 0, 0, 0)