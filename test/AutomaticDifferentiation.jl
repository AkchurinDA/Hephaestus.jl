@testset "Cantilever beam subjected to concentrated load: Linear elastic analysis" begin
    # Define the problem parameters:
    L = 120   # in.
    I = 100   # in.⁴
    E = 29000 # ksi
    P = 100   # kip
    t = float.([L, I, E, P])

    # Define the function to be differentiated:
    function f(t::AbstractVector{<:Real})::Real
        # Extract the problem parameters:
        L, I, E, P = t

        # Define an empty model:
        model = Model()

        # Define the nodes and DOF supports:
        for (i, x) in enumerate(range(0, L, 11))
            if i == 1
                node!(model, i, x, 0, 0, u_x = true, u_y = true, θ_z = true)
            else
                node!(model, i, x, 0, 0)
            end
        end

        # Define the sections:
        section!(model, 1, 1, I, 1, 1)

        # Define the materials:
        material!(model, 1, E, 0.3, 1)

        # Define the elements:
        for i in 1:10
            element!(model, i, i, i + 1, 1, 1)
        end

        # Define the loads:
        concload!(model, 11, 0, -P, 0, 0, 0, 0)

        # Solve the model using a linear elastic analysis:
        solution = solve(model, LinearElasticAnalysis())

        # Extract the vertical displacement of free end of the cantilever beam:
        Δ = getnodaldisplacements(model, solution, 11)[2]

        # Return the vertical displacement of the free end of the cantilever beam:
        return Δ
    end

    # Define the exact solution:
    function f_exact(t::AbstractVector{<:Real})::Real
        # Extract the problem parameters:
        L, I, E, P = t

        # Compute the vertical displacement of the free end of the cantilever beam:
        Δ = -(P * L ^ 3) / (3 * E * I)

        # Return the vertical displacement of the free end of the cantilever beam:
        return Δ
    end

    # Check:
    @test f(t) ≈ f_exact(t) rtol = 1E-3

    # Preallocate the gradient vector:
    ∇f       = similar(t)
    ∇f_exact = similar(t)

    # Compute the gradient vector using ForwardDiff.jl:
    gradient!(f, ∇f, AutoForwardDiff(), t)
    gradient!(f_exact, ∇f_exact, AutoForwardDiff(), t)
    @test ∇f ≈ ∇f_exact rtol = 1E-3

    # Compute the gradient vector using ReverseDiff.jl:
    gradient!(f, ∇f, AutoReverseDiff(), t)
    gradient!(f_exact, ∇f_exact, AutoReverseDiff(), t)
    @test ∇f ≈ ∇f_exact rtol = 1E-3
end

@testset "Cantilever beam subjected to distributed load: Linear elastic analysis" begin
    # Define the problem parameters:
    L = 120   # in.
    I = 100   # in.⁴
    E = 29000 # ksi
    w = 10    # kip / in.
    t = float.([L, I, E, w])

    # Define the function to be differentiated:
    function f(t::AbstractVector{<:Real})::Real
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
        Δ = getnodaldisplacements(model, solution, 11)[2]

        # Return the vertical displacement of the free end of the cantilever beam:
        return Δ
    end

    # Define the exact solution:
    function f_exact(t::AbstractVector{<:Real})::Real
        # Extract the problem parameters:
        L, I, E, w = t

        # Compute the vertical displacement of the free end of the cantilever beam:
        Δ = -(w * L ^ 4) / (8 * E * I)

        # Return the vertical displacement of the free end of the cantilever beam:
        return Δ
    end

    # Check:
    @test f(t) ≈ f_exact(t) rtol = 1E-3

    # Preallocate the gradient vector:
    ∇f       = similar(t)
    ∇f_exact = similar(t)

    # Compute the gradient vector using ForwardDiff.jl:
    gradient!(f, ∇f, AutoForwardDiff(), t)
    gradient!(f_exact, ∇f_exact, AutoForwardDiff(), t)
    @test ∇f ≈ ∇f_exact rtol = 1E-3

    # Compute the gradient vector using ReverseDiff.jl:
    gradient!(f, ∇f, AutoReverseDiff(), t)
    gradient!(f_exact, ∇f_exact, AutoReverseDiff(), t)
    @test ∇f ≈ ∇f_exact rtol = 1E-3
end