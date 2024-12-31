using Hephaestus

# Define the problem parameters:
L = 180   # in.
A = 9.12  # in.²
I = 110   # in.⁴
E = 29000 # ksi
P = 50    # kip
H = 1     # kip
t = float.([L, A, I, E, P, H])

function f(t::AbstractVector{<:Real})::Real
    # Extract the problem parameters:
    L, A, I, E, P, H = t

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
    section!(model, 1, A, I, 1, 1)

    # Define the materials:
    material!(model, 1, E, 0.3, 0)

    # Define the elements:
    for i in 1:10
        element!(model, i, i, i + 1, 1, 1)
    end

    # Define the loads:
    concload!(model, 11, -P, -H, 0, 0, 0, 0)

    solution = solve(model, NonlinearElasticAnalysis(LCM(1 / 100), 100, 100, :standard, 1E-12));

    Δ = getnodaldisplacements(model, solution, 11)[2]

    return Δ
end