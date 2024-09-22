using Test, Hephaestus

@testset "1st-Order Elastic Analysis" begin
    @testset "Cantilever Beam - Concentrated Loading" begin
        M = Model()

        add_node!(M, 1 , 0  * 12, 0, 0)
        add_node!(M, 2 , 1  * 12, 0, 0)
        add_node!(M, 3 , 2  * 12, 0, 0)
        add_node!(M, 4 , 3  * 12, 0, 0)
        add_node!(M, 5 , 4  * 12, 0, 0)
        add_node!(M, 6 , 5  * 12, 0, 0)
        add_node!(M, 7 , 6  * 12, 0, 0)
        add_node!(M, 8 , 7  * 12, 0, 0)
        add_node!(M, 9 , 8  * 12, 0, 0)
        add_node!(M, 10, 9  * 12, 0, 0)
        add_node!(M, 11, 10 * 12, 0, 0)

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

        add_concentrated_load!(M, 11, 0, -1000, 0, 0, 0, 0)

        Solution = solve(M, O1EAnalysis())

        # Test reactions at the fixed end:
        @test Solution.R[1] ≈ 0
        @test Solution.R[2] ≈ 1000
        @test Solution.R[3] ≈ 0
        @test Solution.R[4] ≈ 0
        @test Solution.R[5] ≈ 0
        @test Solution.R[6] ≈ 1000 * 120

        # Test displacements at the free end:
        @test Solution.U[61] ≈ 0
        @test Solution.U[62] ≈ (-1000 * 120 ^ 3) / (3 * 29000 * 1000)
        @test Solution.U[63] ≈ 0
        @test Solution.U[64] ≈ 0
        @test Solution.U[65] ≈ 0
        @test Solution.U[66] ≈ (-1000 * 120 ^ 2) / (2 * 29000 * 1000)
    end

    @testset "Cantilever Beam - Uniform Loading" begin
        M = Model()

        add_node!(M, 1 , 0  * 12, 0, 0)
        add_node!(M, 2 , 1  * 12, 0, 0)
        add_node!(M, 3 , 2  * 12, 0, 0)
        add_node!(M, 4 , 3  * 12, 0, 0)
        add_node!(M, 5 , 4  * 12, 0, 0)
        add_node!(M, 6 , 5  * 12, 0, 0)
        add_node!(M, 7 , 6  * 12, 0, 0)
        add_node!(M, 8 , 7  * 12, 0, 0)
        add_node!(M, 9 , 8  * 12, 0, 0)
        add_node!(M, 10, 9  * 12, 0, 0)
        add_node!(M, 11, 10 * 12, 0, 0)

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

        add_distributed_load!(M, 1 , 0, -100, 0)
        add_distributed_load!(M, 2 , 0, -100, 0)
        add_distributed_load!(M, 3 , 0, -100, 0)
        add_distributed_load!(M, 4 , 0, -100, 0)
        add_distributed_load!(M, 5 , 0, -100, 0)
        add_distributed_load!(M, 6 , 0, -100, 0)
        add_distributed_load!(M, 7 , 0, -100, 0)
        add_distributed_load!(M, 8 , 0, -100, 0)
        add_distributed_load!(M, 9 , 0, -100, 0)
        add_distributed_load!(M, 10, 0, -100, 0)

        Solution = solve(M, O1EAnalysis())

        # Test reactions at the fixed end:
        @test Solution.R[1] ≈ 0
        @test Solution.R[2] ≈ 100 * 120
        @test Solution.R[3] ≈ 0
        @test Solution.R[4] ≈ 0
        @test Solution.R[5] ≈ 0
        @test Solution.R[6] ≈ 100 * 120 ^ 2 / 2

        # Test displacements at the free end:
        @test Solution.U[61] ≈ 0
        @test Solution.U[62] ≈ (-100 * 120 ^ 4) / (8 * 29000 * 1000)
        @test Solution.U[63] ≈ 0
        @test Solution.U[64] ≈ 0
        @test Solution.U[65] ≈ 0
        @test Solution.U[66] ≈ (-100 * 120 ^ 3) / (6 * 29000 * 1000)
    end
end