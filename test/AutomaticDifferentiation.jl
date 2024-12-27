@testset "Cantilever beam subjected to concentrated load: Linear elastic analysis" begin
    # Define the problem parameters:
    L = 120   # in.
    I = 100   # in.⁴
    E = 29000 # ksi
    P = 100   # kip
    t = float.([L, I, E, P])

    # Define the function to be differentiated:
    function f(t)
        # Extract the problem parameters:
        L, I, E, P = t

        # Define an empty model:
        model = Model()

        # Define the nodes and DOF supports:
        for (ID, x) in enumerate(range(0, L, 11))
            if ID == 1
                node!(model, ID, x, 0, 0, u_x = true, u_y = true, θ_z = true)
            else
                node!(model, ID, x, 0, 0)
            end
        end

        # Define the sections:
        section!(model, 1, 1, I, 1, 1)

        # Define the materials:
        material!(model, 1, E, 0.3, 1)

        # Define the elements:
        for ID in 1:10
            element!(model, ID, ID, ID + 1, 1, 1)
        end

        # Define the loads:
        concload!(model, 11, 0, -P, 0, 0, 0, 0)

        # Solve the model using a linear elastic analysis:
        solution = solve(model, LinearElasticAnalysis())

        # Extract the vertical displacement of free end of the cantilever beam:
        Δ = -getnodaldisplacements(model, solution, 11)[2]

        # Return the vertical displacement of the free end of the cantilever beam:
        return Δ
    end

    # Check:
    @test f(t) ≈ (P * L ^ 3) / (3 * E * I) rtol = 1E-3

    # Define the exact gradient vector:
    ∇f = [
        +(P * L ^ 2) / (E * I),
        -((P * L ^ 3) / (3 * E * I ^ 2)),
        -((P * L ^ 3) / (3 * E ^ 2 * I)),
        +(L ^ 3) / (3 * E * I)]

    # Preallocate the gradient vector:
    ∇f′ = similar(t)

    # Compute the gradient vector using ForwardDiff.jl:
    gradient!(f, ∇f′, AutoForwardDiff(), t)
    @test ∇f ≈ ∇f′ rtol = 1E-3

    # Compute the gradient vector using ReverseDiff.jl:
    gradient!(f, ∇f′, AutoReverseDiff(), t)
    @test ∇f ≈ ∇f′ rtol = 1E-3
end

@testset "Cantilever beam subjected to distributed load: Linear elastic analysis" begin
    # Define the problem parameters:
    L = 120   # in.
    I = 100   # in.⁴
    E = 29000 # ksi
    w = 10    # kip / in.
    t = float.([L, I, E, w])

    # Define the function to be differentiated:
    function f(t)
        # Extract the problem parameters:
        L, I, E, w = t

        # Define an empty model:
        model = Model()

        # Define the nodes and DOF supports:
        for (ID, x) in enumerate(range(0, L, 11))
            if ID == 1
                node!(model, ID, x, 0, 0, u_x = true, u_y = true, θ_z = true)
            else
                node!(model, ID, x, 0, 0)
            end
        end

        # Define the sections:
        section!(model, 1, 1, I, 1, 1)

        # Define the materials:
        material!(model, 1, E, 0.3, 1)

        # Define the elements:
        for ID in 1:10
            element!(model, ID, ID, ID + 1, 1, 1)
        end

        # Define the loads:
        for ID in 1:10
            distload!(model, ID, 0, -w, 0)
        end

        # Solve the model using a linear elastic analysis:
        solution = solve(model, LinearElasticAnalysis())

        # Extract the vertical displacement of free end of the cantilever beam:
        Δ = -getnodaldisplacements(model, solution, 11)[2]

        # Return the vertical displacement of the free end of the cantilever beam:
        return Δ
    end

    # Check:
    @test f(t) ≈ (w * L ^ 4) / (8 * E * I) rtol = 1E-3

    # Define the exact gradient vector:
    ∇f = [
        +(w * L ^ 3) / (2 * E * I),
        -((w * L ^ 4) / (8 * E * I ^ 2)),
        -((w * L ^ 4) / (8 * E ^ 2 * I)),
        +(L ^ 4) / (8 * E * I)]

    # Preallocate the gradient vector:
    ∇f′ = similar(t)

    # Compute the gradient vector using ForwardDiff.jl:
    gradient!(f, ∇f′, AutoForwardDiff(), t)
    @test ∇f ≈ ∇f′ rtol = 1E-3

    # Compute the gradient vector using ReverseDiff.jl:
    gradient!(f, ∇f′, AutoReverseDiff(), t)
    @test ∇f ≈ ∇f′ rtol = 1E-3
end