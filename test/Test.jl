using Hephaestus
using ForwardDiff

function f(x)
    # Create a new model:
    M = Model()

    # Add nodes:
    add_node!(M, 1, (1 - 1) * 1000, 0, 0)
    add_node!(M, 2, (2 - 1) * 1000, 0, 0)
    add_node!(M, 3, (3 - 1) * 1000, 0, 0)
    add_node!(M, 4, (4 - 1) * 1000, 0, 0)
    add_node!(M, 5, (5 - 1) * 1000, 0, 0)
    add_node!(M, 6, (6 - 1) * 1000, 0, 0)

    # Add materials:
    add_material!(M, 1, 200_000, 0.3, 9.81 * 7850)

    # Add sections:
    add_section!(M, 1, 4E3, 50E6, 0, 0)

    # Add elements:
    add_element!(M, 1, 1, 2, 1, 1)
    add_element!(M, 2, 2, 3, 1, 1)
    add_element!(M, 3, 3, 4, 1, 1)
    add_element!(M, 4, 4, 5, 1, 1)
    add_element!(M, 5, 5, 6, 1, 1)

    add_nodal_load!(M, 6, 0, x, 0, 0, 0, 0)

    # Add supports:
    add_support!(M, 1, true, true , true, true, true, true )
    add_support!(M, 2, true, false, true, true, true, false)
    add_support!(M, 3, true, false, true, true, true, false)
    add_support!(M, 4, true, false, true, true, true, false)
    add_support!(M, 5, true, false, true, true, true, false)
    add_support!(M, 6, true, false, true, true, true, false)

    U = solve(M, :O1E)

    return U[end]
end

@time ForwardDiff.derivative(f, -1000)