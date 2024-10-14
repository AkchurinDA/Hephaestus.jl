using Test
using Hephaestus
# using DifferentiationInterface, ForwardDiff, ReverseDiff

@testset "Structural Analysis: 1st-Order Elastic Analysis" begin
    @testset "Cantilever Beam - Concentrated Loading" begin
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

        add_material!(model, 1, 29000, 0.3, 0.000284)

        add_section!(model, 1, 5, 1000, 1000, 10)

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

        add_support!(model, 1 , true , true , true, true, true, true )
        add_support!(model, 2 , false, false, true, true, true, false)
        add_support!(model, 3 , false, false, true, true, true, false)
        add_support!(model, 4 , false, false, true, true, true, false)
        add_support!(model, 5 , false, false, true, true, true, false)
        add_support!(model, 6 , false, false, true, true, true, false)
        add_support!(model, 7 , false, false, true, true, true, false)
        add_support!(model, 8 , false, false, true, true, true, false)
        add_support!(model, 9 , false, false, true, true, true, false)
        add_support!(model, 10, false, false, true, true, true, false)
        add_support!(model, 11, false, false, true, true, true, false)

        add_conc_load!(model, 11, 0, -1000, 0, 0, 0, 0)

        solution = solve(model, LinearElasticAnalysis())

        # Test reactions at the fixed end:
        @test solution.R[1] ≈ 0
        @test solution.R[2] ≈ 1000
        @test solution.R[3] ≈ 0
        @test solution.R[4] ≈ 0
        @test solution.R[5] ≈ 0
        @test solution.R[6] ≈ 1000 * 120

        # Test displacements at the free end:
        @test solution.U[61] ≈ 0
        @test solution.U[62] ≈ (-1000 * 120 ^ 3) / (3 * 29000 * 1000)
        @test solution.U[63] ≈ 0
        @test solution.U[64] ≈ 0
        @test solution.U[65] ≈ 0
        @test solution.U[66] ≈ (-1000 * 120 ^ 2) / (2 * 29000 * 1000)
    end

    @testset "Cantilever Beam - Uniform Loading" begin
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

        add_material!(model, 1, 29000, 0.3, 0.000284)

        add_section!(model, 1, 5, 1000, 1000, 10)

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

        add_support!(model, 1 , true , true , true, true, true, true )
        add_support!(model, 2 , false, false, true, true, true, false)
        add_support!(model, 3 , false, false, true, true, true, false)
        add_support!(model, 4 , false, false, true, true, true, false)
        add_support!(model, 5 , false, false, true, true, true, false)
        add_support!(model, 6 , false, false, true, true, true, false)
        add_support!(model, 7 , false, false, true, true, true, false)
        add_support!(model, 8 , false, false, true, true, true, false)
        add_support!(model, 9 , false, false, true, true, true, false)
        add_support!(model, 10, false, false, true, true, true, false)
        add_support!(model, 11, false, false, true, true, true, false)

        add_dist_load!(model, 1 , 0, -100, 0)
        add_dist_load!(model, 2 , 0, -100, 0)
        add_dist_load!(model, 3 , 0, -100, 0)
        add_dist_load!(model, 4 , 0, -100, 0)
        add_dist_load!(model, 5 , 0, -100, 0)
        add_dist_load!(model, 6 , 0, -100, 0)
        add_dist_load!(model, 7 , 0, -100, 0)
        add_dist_load!(model, 8 , 0, -100, 0)
        add_dist_load!(model, 9 , 0, -100, 0)
        add_dist_load!(model, 10, 0, -100, 0)

        solution = solve(model, LinearElasticAnalysis())

        # Test reactions at the fixed end:
        @test solution.R[1] ≈ 0
        @test solution.R[2] ≈ 100 * 120
        @test solution.R[3] ≈ 0
        @test solution.R[4] ≈ 0
        @test solution.R[5] ≈ 0
        @test solution.R[6] ≈ 100 * 120 ^ 2 / 2

        # Test displacements at the free end:
        @test solution.U[61] ≈ 0
        @test solution.U[62] ≈ (-100 * 120 ^ 4) / (8 * 29000 * 1000)
        @test solution.U[63] ≈ 0
        @test solution.U[64] ≈ 0
        @test solution.U[65] ≈ 0
        @test solution.U[66] ≈ (-100 * 120 ^ 3) / (6 * 29000 * 1000)
    end
end