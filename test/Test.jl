using Hephaestus
# using FiniteDiff
# using ForwardDiff
# using BenchmarkTools

M = Model()

add_node!(M, 1, 0 * 12, 0, 0)
add_node!(M, 2, 1 * 12, 0, 0)
add_node!(M, 3, 2 * 12, 0, 0)
add_node!(M, 4, 3 * 12, 0, 0)
add_node!(M, 5, 4 * 12, 0, 0)
add_node!(M, 6, 5 * 12, 0, 0)

add_material!(M, 1, 29000, 0.3, 0.000284)

add_section!(M, 1, 5, 1000, 1000, 10)

add_element!(M, 1, 1, 2, 1, 1, 0)
add_element!(M, 2, 2, 3, 1, 1, 0)
add_element!(M, 3, 3, 4, 1, 1, 0)
add_element!(M, 4, 4, 5, 1, 1, 0)
add_element!(M, 5, 5, 6, 1, 1, 0)

add_support!(M, 1, true , true, true, true, true, true )
add_support!(M, 2, false, true, true, true, true, false)
add_support!(M, 3, false, true, true, true, true, false)
add_support!(M, 4, false, true, true, true, true, false)
add_support!(M, 5, false, true, true, true, true, false)
add_support!(M, 6, false, true, true, true, true, false)

Solution = solve(M, EB())

# @benchmark FiniteDiff.finite_difference_derivative(f, -1000.0)
# @benchmark ForwardDiff.derivative(f, -1000.0)