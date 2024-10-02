using Hephaestus
using GLMakie

M = Model()

add_node!(M, 1 , 0, 0  * 12, 0)
add_node!(M, 2 , 0, 1  * 12, 0)
add_node!(M, 3 , 0, 2  * 12, 0)
add_node!(M, 4 , 0, 3  * 12, 0)
add_node!(M, 5 , 0, 4  * 12, 0)
add_node!(M, 6 , 0, 5  * 12, 0)
add_node!(M, 7 , 0, 6  * 12, 0)
add_node!(M, 8 , 0, 7  * 12, 0)
add_node!(M, 9 , 0, 8  * 12, 0)
add_node!(M, 10, 0, 9  * 12, 0)
add_node!(M, 11, 0, 10 * 12, 0)

add_material!(M, 1, 29000, 0.3, 0.000284)

add_section!(M, 1, 5, 1000, 1000, 10)

add_element!(M, 1 , 1 , 2 , 1, 1)
add_element!(M, 2 , 2 , 3 , 1, 1)
add_element!(M, 3 , 3 , 4 , 1, 1)
add_element!(M, 4 , 4 , 5 , 1, 1)
add_element!(M, 5 , 5 , 6 , 1, 1)
add_element!(M, 6 , 6 , 7 , 1, 1)
add_element!(M, 7 , 7 , 8 , 1, 1)
add_element!(M, 8 , 8 , 9 , 1, 1)
add_element!(M, 9 , 9 , 10, 1, 1)
add_element!(M, 10, 10, 11, 1, 1)

add_support!(M, 1 , true , true , true, true, true, true )
add_support!(M, 2 , false, false, true, true, true, false)
add_support!(M, 3 , false, false, true, true, true, false)
add_support!(M, 4 , false, false, true, true, true, false)
add_support!(M, 5 , false, false, true, true, true, false)
add_support!(M, 6 , false, false, true, true, true, false)
add_support!(M, 7 , false, false, true, true, true, false)
add_support!(M, 8 , false, false, true, true, true, false)
add_support!(M, 9 , false, false, true, true, true, false)
add_support!(M, 10, false, false, true, true, true, false)
add_support!(M, 11, false, false, true, true, true, false)

add_conc_load!(M, 11, 0, -1, 0, 0, 0, 0)

plotmodel(M)

begin
    F = Figure()
    A = Axis3(F[1, 1])
    
    plotmodel!(A, M)

    display(F)
end