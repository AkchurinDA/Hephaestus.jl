@testset "Model building" begin
    # Define an empty model:
    model = Model()

    # Define the nodes and DOF supports: 
    node!(model, 1,   0, 0, 0)
    node!(model, 2, 120, 0, 0)

    # Test the node function:
    @test_throws AssertionError node!(model, 1,  0, 0, 0) # Duplicate node (i)
    @test_throws AssertionError node!(model, 2, 12, 0, 0) # Duplicate node (j)

    # Define the sections:
    section!(model, 1, 10, 100, 100, 5)

    # Test the section function:
    @test_throws AssertionError section!(model, 1, 10, 100, 100, 5) # Duplicate section

    # Define the materials:
    material!(model, 1, 29000, 0.3, 0.290)

    # Test the material function:
    @test_throws AssertionError material!(model, 1, 29000, 0.3, 0.290) # Duplicate material

    # Define the elements:
    element!(model, 1, 1, 2, 1, 1)

    # Test the element function:
    @test_throws AssertionError element!(model, 1, 1, 2, 1, 1) # Duplicate element
    @test_throws AssertionError element!(model, 2, 1, 4, 1, 1) # Invalid node (i)
    @test_throws AssertionError element!(model, 2, 3, 2, 1, 1) # Invalid node (j)
    @test_throws AssertionError element!(model, 2, 3, 4, 2, 1) # Invalid section
    @test_throws AssertionError element!(model, 2, 3, 4, 1, 2) # Invalid material

    # Define the concentrated loads:
    concload!(model, 1, 1, 0, 0, 0, 0, 0)

    # Test the concentrated load function:
    @test_throws AssertionError concload!(model, 1, 1, 0, 0, 0, 0, 0) # Duplicate concentrated load
    @test_throws AssertionError concload!(model, 3, 1, 0, 0, 0, 0, 0) # Invalid node
    @test_throws AssertionError concload!(model, 2, 0, 0, 0, 0, 0, 0) # Invalid load

    # Define the distributed loads:
    distload!(model, 1, 0, -1, 0)

    # Test the distributed load function:
    @test_throws AssertionError distload!(model, 1, 0, -1, 0) # Duplicate distributed load
    @test_throws AssertionError distload!(model, 2, 0, -1, 0) # Invalid element
end