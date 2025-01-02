abstract type AbstractNonlinearSolver end

"""
    struct LCM

A type representing the load control method.
"""
struct LCM <: AbstractNonlinearSolver
    # Load increment:
    Δλ::Real
end

function coefficients(nonlinearsolver::LCM, partitionindices::Vector{Bool})
    a = zeros(length(partitionindices))
    b = 1
    c = nonlinearsolver.Δλ

    return a, b, c
end

"""
    struct DCM

A type representing the displacement control method.
"""
struct DCM <: AbstractNonlinearSolver
    # Displacement increment:
    Δu::Real
end

"""
    struct WCM

A type representing the work control method.
"""
struct WCM <: AbstractNonlinearSolver
    # Work increment:
    ΔW::Real
end

"""
    struct ALCM

A type representing the arc-length control method.
"""
struct ALCM <: AbstractNonlinearSolver
    # Arc-length increment:
    ΔS::Real
end

function constraintequation(
    a   ::AbstractVector{<:Real},
    b   ::Real,
    c   ::Real,
    δu_p::AbstractVector{<:Real},
    δu_r::AbstractVector{<:Real})::Real
    # Compute the load paramter increment:
    δλ = (c - dot(a, δu_r)) / (b + dot(a, δu_p))

    # Return the load parameter increment:
    return δλ
end
