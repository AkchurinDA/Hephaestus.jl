"""
    struct Element

A type that represents an element in the FE model of a structure.

This type should not be used directly. 
Instead, use the [`add_element!()`](@ref) function to add an element to the model. 
To remove an element from the model, use the [`del_element!()`](@ref) function.
"""
mutable struct Element{CTI <: Real, CTJ <: Real, MPT <: Real, SPT <: Real}
    "Unique identifier of the element provided by the user"
    ID_u        ::Int
    "Unique identifier of the element assigned by the program"
    ID_i        ::Union{Nothing, Int}
    "ID of the 1st node"
    node_i_ID   ::Int
    "``x``-coordinate of the 1st node"
    x_i         ::CTI
    "``y``-coordinate of the 1st node"
    y_i         ::CTI
    "``z``-coordinate of the 1st node"
    z_i         ::CTI
    "ID of the 2nd node"
    node_j_ID   ::Int
    "``x``-coordinate of the 2nd node"
    x_j         ::CTJ
    "``y``-coordinate of the 2nd node"
    y_j         ::CTJ
    "``z``-coordinate of the 2nd node"
    z_j         ::CTJ
    "ID of the material"
    material_ID ::Int
    "Young's modulus"
    E           ::MPT
    "Poisson's ratio"
    ν           ::MPT
    "Density"
    ρ           ::MPT
    "ID of the section"
    section_ID  ::Int
    "Cross-sectional area of the element"
    A           ::SPT
    "Moment of inertia about the local z-axis"
    I_zz        ::SPT
    "Moment of inertia about the local y-axis"
    I_yy        ::SPT
    "Polar moment of inertia"
    J           ::SPT
    "Angle that defines the orientation of section's local y-axis"
    ω           ::Real
    "DOF releases at the ends of the element"
    releases    ::Vector{Bool}
    "Length of the element"
    L           ::Real
    "Transformation matrix"
    T           ::Matrix{<:Real}
    "Uncondensed element elastic stiffness matrix in its local coordinate system"
    k_e_l       ::Matrix{<:Real}
    "Uncondensed element elastic stiffness matrix in the global coordinate system"
    k_e_g       ::Matrix{<:Real}
    "Condensed element elastic stiffness matrix in its local coordinate system"
    k_e_l_c     ::Matrix{<:Real}
    "Condensed element elastic stiffness matrix in the global coordinate system"
    k_e_g_c     ::Matrix{<:Real}
end

function _compute_L(
    x_i::Real, y_i::Real, z_i::Real,
    x_j::Real, y_j::Real, z_j::Real)
    L = sqrt((x_j - x_i) ^ 2 + (y_j - y_i) ^ 2 + (z_j - z_i) ^ 2)

    return L
end

function _compute_T(
    # Nodal coordinates of the element's 1st node:
    x_i::CTI, y_i::CTI, z_i::CTI,
    # Nodal coordinates of the element's 2nd node:
    x_j::CTJ, y_j::CTJ, z_j::CTJ,
    # Element properties:
    L::EPT,
    ω::Real) where {CTI <: Real, CTJ <: Real, EPT <: Real}
    # Preallocate:
    S = float(promote_type(CTI, CTJ, EPT, typeof(ω)))
    T = zeros(S, 12, 12)

    # Construct the sub-rotation matrix:
    ρ   = -atan(z_j - z_i, x_j - x_i)
    χ   = π / 2 - acos((y_j - y_i) / L)
    γ_ρ = [
        +cos(ρ)  0 -sin(ρ); 
        0       +1       0; 
        +sin(ρ)  0 +cos(ρ)]
    γ_χ = [
        +cos(χ) +sin(χ)  0; 
        -sin(χ) +cos(χ)  0; 
            0       0   +1]
    γ_ω = [     
        +1       0       0;
         0 +cos(ω) +sin(ω); 
         0 -sin(ω) +cos(ω)]
    γ = γ_ω * γ_χ * γ_ρ

    # Remove small values if any:
    γ[abs.(γ) .< eps()] .= 0

    # Define the components of the transformation matrix:
    @inbounds T[1:3  , 1:3  ] .= γ
    @inbounds T[4:6  , 4:6  ] .= γ
    @inbounds T[7:9  , 7:9  ] .= γ
    @inbounds T[10:12, 10:12] .= γ

    # Return the transformation matrix:
    return T
end

function _compute_k_e_l(
    # Element properties:
    L::EPT,
    # Material properties:
    E::MPT, ν::MPT,
    # Section properties:
    A::SPT, I_zz::SPT, I_yy::SPT, J::SPT) where {EPT <: Real, MPT <: Real, SPT <: Real}
    # Preallocate:
    T     = float(promote_type(EPT, MPT, SPT))
    k_e_l = zeros(T, 12, 12)

    # Compute the components of the element elastic stiffness matrix in its upper triangular part:
    @inbounds k_e_l[1 , 1 ] = +E * A / L
    @inbounds k_e_l[1 , 7 ] = -E * A / L
    @inbounds k_e_l[2 , 2 ] = +12 * E * I_zz / L^3
    @inbounds k_e_l[2 , 6 ] = +6 * E * I_zz / L^2
    @inbounds k_e_l[2 , 8 ] = -12 * E * I_zz / L^3
    @inbounds k_e_l[2 , 12] = +6 * E * I_zz / L^2
    @inbounds k_e_l[3 , 3 ] = +12 * E * I_yy / L^3
    @inbounds k_e_l[3 , 5 ] = -6 * E * I_yy / L^2
    @inbounds k_e_l[3 , 9 ] = -12 * E * I_yy / L^3
    @inbounds k_e_l[3 , 11] = -6 * E * I_yy / L^2
    @inbounds k_e_l[4 , 4 ] = +E * J / (2 * (1 + ν) * L)
    @inbounds k_e_l[4 , 10] = -E * J / (2 * (1 + ν) * L)
    @inbounds k_e_l[5 , 5 ] = +4 * E * I_yy / L
    @inbounds k_e_l[5 , 9 ] = +6 * E * I_yy / L^2
    @inbounds k_e_l[5 , 11] = +2 * E * I_yy / L
    @inbounds k_e_l[6 , 6 ] = +4 * E * I_zz / L
    @inbounds k_e_l[6 , 8 ] = -6 * E * I_zz / L^2
    @inbounds k_e_l[6 , 12] = +2 * E * I_zz / L
    @inbounds k_e_l[7 , 7 ] = +E * A / L
    @inbounds k_e_l[8 , 8 ] = +12 * E * I_zz / L^3
    @inbounds k_e_l[8 , 12] = -6 * E * I_zz / L^2
    @inbounds k_e_l[9 , 9 ] = +12 * E * I_yy / L^3
    @inbounds k_e_l[9 , 11] = +6 * E * I_yy / L^2
    @inbounds k_e_l[10, 10] = +E * J / (2 * (1 + ν) * L)
    @inbounds k_e_l[11, 11] = +4 * E * I_yy / L
    @inbounds k_e_l[12, 12] = +4 * E * I_zz / L

    # Compute the components of the element elastic stiffness matrix in its lower triangular part:
    for i in 1:12, j in (i + 1):12
        @inbounds k_e_l[j, i] = k_e_l[i, j]
    end

    # Remove small values if any:
    k_e_l[abs.(k_e_l) .< eps()] .= 0

    # Return the element elastic stiffness matrix:
    return k_e_l
end

function _compute_k_g_l(
    # Element properties:
    L::EPT,
    # Section properties:
    A::SPT, I_zz::SPT, I_yy::SPT,
    # Axial load:
    P::ALT) where {EPT <: Real, SPT <: Real, ALT <: Real}
    # Preallocate:
    T         = float(promote_type(EPT, SPT, ALT))
    k_g_l = zeros(T, 12, 12)

    # Compute the components of the element geometric stiffness matrix in its upper triangular part:
    @inbounds k_g_l[1 , 1 ] = +P / L
    @inbounds k_g_l[1 , 7 ] = -P / L
    @inbounds k_g_l[2 , 2 ] = +6 * P / (5 * L)
    @inbounds k_g_l[2 , 6 ] = +P / 10
    @inbounds k_g_l[2 , 8 ] = -6 * P / (5 * L)
    @inbounds k_g_l[2 , 12] = +P / 10
    @inbounds k_g_l[3 , 3 ] = +6 * P / (5 * L)
    @inbounds k_g_l[3 , 5 ] = -P / 10
    @inbounds k_g_l[3 , 9 ] = -6 * P / (5 * L)
    @inbounds k_g_l[3 , 11] = -P / 10
    @inbounds k_g_l[4 , 4 ] = +P * (I_zz + I_yy) / (A * L)
    @inbounds k_g_l[4 , 10] = -P * (I_zz + I_yy) / (A * L)
    @inbounds k_g_l[5 , 5 ] = +2 * P * L / 15
    @inbounds k_g_l[5 , 9 ] = +P / 10
    @inbounds k_g_l[5 , 11] = -P * L / 30
    @inbounds k_g_l[6 , 6 ] = +2 * P * L / 15
    @inbounds k_g_l[6 , 8 ] = -P / 10
    @inbounds k_g_l[6 , 12] = -P * L / 30
    @inbounds k_g_l[7 , 7 ] = +P / L
    @inbounds k_g_l[8 , 8 ] = +6 * P / (5 * L)
    @inbounds k_g_l[8 , 12] = -P / 10
    @inbounds k_g_l[9 , 9 ] = +6 * P / (5 * L)
    @inbounds k_g_l[9 , 11] = +P / 10
    @inbounds k_g_l[10, 10] = +P * (I_zz + I_yy) / (A * L)
    @inbounds k_g_l[11, 11] = +2 * P * L / 15
    @inbounds k_g_l[12, 12] = +2 * P * L / 15

    # Compute the components of the element geometric stiffness matrix in its lower triangular part:
    for i in 1:12, j in (i + 1):12
        @inbounds k_g_l[j, i] = k_g_l[i, j]
    end

    # Remove small values if any:
    k_g_l[abs.(k_g_l) .< eps()] .= 0

    # Return the element geometric stiffness matrix:
    return k_g_l
end

function _condense(k::Matrix{<:Real}, releases::Vector{Bool})
    # Find released DOFs:
    released_DOFs = findall(x -> x == true, releases) # Released DOFs

    # Copy the stiffness matrix:
    k_c = copy(k)

    # Set rows and columns corresponding to released DOFs to zero:
    k_c[released_DOFs, :] .= 0 # Set rows corresponding to released DOFs to zero
    k_c[:, released_DOFs] .= 0 # Set columns corresponding to released DOFs to zero

    # Return the condensed stiffness matrix:
    return k_c
end