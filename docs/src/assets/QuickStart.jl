using Hephaestus
using CairoMakie
CairoMakie.activate!(type = "svg")
CairoMakie.set_theme!(theme_latexfonts())

model = Model()

add_node!(model, 1 , 0  * 12, 0, 0)
add_node!(model, 2 , 1  * 12, 0, 0)
add_node!(model, 3 , 2  * 12, 0, 0)
add_node!(model, 4 , 3  * 12, 0, 0)
add_node!(model, 5 , 4  * 12, 0, 0)
add_node!(model, 6 , 5  * 12, 0, 0)
add_node!(model, 7 , 6  * 12, 0, 0)
add_node!(model, 8 , 7  * 12, 0, 0)
add_node!(model, 9 , 8  * 12, 0, 0)
add_node!(model, 10, 9  * 12, 0, 0)
add_node!(model, 11, 10 * 12, 0, 0)

add_material!(model, 1, 29000, 0.3, 0.000290)

add_section!(model, 1, 9.16, 180, 180, 359)

add_element!(model, 1 , 1 , 2 , 1, 1)
add_element!(model, 2 , 2 , 3 , 1, 1)
add_element!(model, 3 , 3 , 4 , 1, 1)
add_element!(model, 4 , 4 , 5 , 1, 1)
add_element!(model, 5 , 5 , 6 , 1, 1)
add_element!(model, 6 , 6 , 7 , 1, 1)
add_element!(model, 7 , 7 , 8 , 1, 1)
add_element!(model, 8 , 8 , 9 , 1, 1)
add_element!(model, 9 , 9 , 10, 1, 1)
add_element!(model, 10, 10, 11, 1, 1)

add_support!(model, 1 , true , true, true , true, true , true)
add_support!(model, 2 , false, true, false, true, false, true)
add_support!(model, 3 , false, true, false, true, false, true)
add_support!(model, 4 , false, true, false, true, false, true)
add_support!(model, 5 , false, true, false, true, false, true)
add_support!(model, 6 , false, true, false, true, false, true)
add_support!(model, 7 , false, true, false, true, false, true)
add_support!(model, 8 , false, true, false, true, false, true)
add_support!(model, 9 , false, true, false, true, false, true)
add_support!(model, 10, false, true, false, true, false, true)
add_support!(model, 11, false, true, false, true, false, true)

add_dist_load!(model, 1 , 0, 0, -1)
add_dist_load!(model, 2 , 0, 0, -1)
add_dist_load!(model, 3 , 0, 0, -1)
add_dist_load!(model, 4 , 0, 0, -1)
add_dist_load!(model, 5 , 0, 0, -1)
add_dist_load!(model, 6 , 0, 0, -1)
add_dist_load!(model, 7 , 0, 0, -1)
add_dist_load!(model, 8 , 0, 0, -1)
add_dist_load!(model, 9 , 0, 0, -1)
add_dist_load!(model, 10, 0, 0, -1)

solution = solve(model, O1EAnalysis())

begin
    F = Figure(size = 72 .* (12, 6))

    A = Axis3(F[1, 1],
        title       = "Model",
        xlabel      = L"$x$ (in.)", xlabelrotation = 0,
        ylabel      = L"$y$ (in.)", ylabelrotation = 0,
        zlabel      = L"$z$ (in.)", zlabelrotation = +π / 2,
        limits      = (nothing, nothing, -60, +60, -60, +60),
        aspect      = :data,
        azimuth     = -π / 4, elevation = +π / 12,
        protrusions = 60)

    plotmodel!(A, model)

    A = Axis3(F[1, 2],
        title       = "Solution",
        xlabel      = L"$x$ (in.)", xlabelrotation = 0,
        ylabel      = L"$y$ (in.)", ylabelrotation = 0,
        zlabel      = L"$z$ (in.)", zlabelrotation = +π / 2,
        limits      = (nothing, nothing, -60, +60, -60, +60),
        aspect      = :data,
        azimuth     = -π / 4, elevation = +π / 12,
        protrusions = 60)

    plotsolution!(A, model, solution,
        scale = 5)

    display(F)
end

save("docs/src/assets/QuickStart.pdf", F)