abstract type AbstractNonlinearSolver end

struct LoadControl <: AbstractNonlinearSolver
    Δλ::Real
end

struct ArcLengthControl <: AbstractNonlinearSolver
    Δs::Real
end

function initialize(model::Model, nonlinearsolver::AbstractNonlinearSolver)
    # a = ?
    # b = ?
    # c = ?
end

function constraint(a, b, c, δu_p, δu_r)
    # Compute the load paramter increment:
    δλ = (c - dot(a, δu_r)) / (b + dot(a, δu_p))

    # Return the load parameter increment:
    return δλ
end