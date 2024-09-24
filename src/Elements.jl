"""
    struct Element

A type representing an element in the finite element in a structure of interest.
"""
struct Element{CTI<:Real, CTJ<:Real, MPT<:Real, SPT<:Real}
    ID            ::Int
    node_i_ID     ::Int
    x_i           ::CTI
    y_i           ::CTI
    z_i           ::CTI
    node_j_ID     ::Int
    x_j           ::CTJ
    y_j           ::CTJ
    z_j           ::CTJ
    E             ::MPT
    ν             ::MPT
    ρ             ::MPT
    A             ::SPT
    I_zz          ::SPT
    I_yy          ::SPT
    J             ::SPT
    ω             ::Real
    releases_i    ::Vector{Bool}
    releases_j    ::Vector{Bool}
    L             ::Real
    γ             ::Matrix{<:Real}
    T             ::BlockDiagonal{<:Real}
    k_e_l         ::Matrix{<:Real}
    k_e_g         ::Matrix{<:Real}
    k_g_l         ::Matrix{<:Real}
    k_g_g         ::Matrix{<:Real}
end

function _compute_L(
    x_i::CTI, y_i::CTI, z_i::CTI,
    x_j::CTJ, y_j::CTJ, z_j::CTJ) where {CTI<:Real, CTJ<:Real}
    # Compute the length of the element:
    L = sqrt((x_j - x_i) ^ 2 + (y_j - y_i) ^ 2 + (z_j - z_i) ^ 2)

    # Return the length of the element:
    return L
end

function _compute_T(
    x_i::CTI, y_i::CTI, z_i::CTI,
    x_j::CTJ, y_j::CTJ, z_j::CTJ,
    L::Real,  
    ω::Real) where {CTI<:Real, CTJ<:Real}
    # Compute the rotation angles:
    ρ = -atan(z_j - z_i, x_j - x_i)
    χ = π / 2 - acos((y_j - y_i) / L)

    # Construct the sub-transformation matrix:
    s_ρ, c_ρ = sincos(ρ)
    s_χ, c_χ = sincos(χ)
    s_ω, c_ω = sincos(ω)
    γ = [
        +c_χ * c_ρ                      +s_χ          -c_χ * s_ρ                  
        +s_ω * s_ρ - c_ω * s_χ * c_ρ    +c_ω * c_χ    +c_ω * s_χ * s_ρ + s_ω * c_ρ
        +c_ω * s_ρ + s_ω * s_χ * c_ρ    -s_ω * c_χ    -s_ω * s_χ * s_ρ + c_ω * c_ρ]

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, γ, γ)

    # Preallocate the transformation matrix and fill it:
    T = BlockDiagonal([γ, γ, γ, γ])

    # Return the transformation matrices:
    return γ, T
end

function _compute_k_e_l(
    E::MPT, ν::MPT,
    A::SPT, I_zz::SPT, I_yy::SPT, J::SPT,
    L::EPT) where {MPT<:Real, SPT<:Real, EPT<:Real}
    # Preallocate:
    T     = float(promote_type(MPT, SPT, EPT))
    k_e_l = zeros(T, 12, 12)

    # Precompute the element stiffnesses:
    k_e_a    = E * A / L
    k_e_b_zz = E * I_zz / L
    k_e_b_yy = E * I_yy / L
    k_e_t    = E * J / (2 * (1 + ν) * L)

    # Compute the components of the element elastic stiffness matrix in its upper triangular part:
    @inbounds k_e_l[1 , 1 ] = +k_e_a
    @inbounds k_e_l[1 , 7 ] = -k_e_a
    @inbounds k_e_l[2 , 2 ] = +12 * k_e_b_zz / L ^ 2
    @inbounds k_e_l[2 , 6 ] = +6  * k_e_b_zz / L
    @inbounds k_e_l[2 , 8 ] = -12 * k_e_b_zz / L ^ 2
    @inbounds k_e_l[2 , 12] = +6  * k_e_b_zz / L
    @inbounds k_e_l[3 , 3 ] = +12 * k_e_b_yy / L ^ 2
    @inbounds k_e_l[3 , 5 ] = -6  * k_e_b_yy / L
    @inbounds k_e_l[3 , 9 ] = -12 * k_e_b_yy / L ^ 2
    @inbounds k_e_l[3 , 11] = -6  * k_e_b_yy / L
    @inbounds k_e_l[4 , 4 ] = +k_e_t
    @inbounds k_e_l[4 , 10] = -k_e_t
    @inbounds k_e_l[5 , 5 ] = +4 * k_e_b_yy
    @inbounds k_e_l[5 , 9 ] = +6 * k_e_b_yy / L
    @inbounds k_e_l[5 , 11] = +2 * k_e_b_yy
    @inbounds k_e_l[6 , 6 ] = +4 * k_e_b_zz
    @inbounds k_e_l[6 , 8 ] = -6 * k_e_b_zz / L
    @inbounds k_e_l[6 , 12] = +2 * k_e_b_zz
    @inbounds k_e_l[7 , 7 ] = +k_e_a
    @inbounds k_e_l[8 , 8 ] = +12 * k_e_b_zz / L ^ 2
    @inbounds k_e_l[8 , 12] = -6  * k_e_b_zz / L
    @inbounds k_e_l[9 , 9 ] = +12 * k_e_b_yy / L ^ 2
    @inbounds k_e_l[9 , 11] = +6  * k_e_b_yy / L
    @inbounds k_e_l[10, 10] = +k_e_t
    @inbounds k_e_l[11, 11] = +4 * k_e_b_yy
    @inbounds k_e_l[12, 12] = +4 * k_e_b_zz

    # Compute the components of the element elastic stiffness matrix in its lower triangular part:
    for i in 1:12, j in (i + 1):12
        @inbounds k_e_l[j, i] = k_e_l[i, j]
    end

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, k_e_l, k_e_l)
    
    # Return the element elastic stiffness matrix:
    return k_e_l
end

function _compute_k_g_l(
    A::SPT, I_zz::SPT, I_yy::SPT,
    L::EPT) where {SPT<:Real, EPT<:Real}
    # Preallocate:
    T     = float(promote_type(SPT, EPT))
    k_g_l = zeros(T, 12, 12)

    # Compute the components of the element geometric stiffness matrix in its upper triangular part:
    @inbounds k_g_l[1 , 1 ] = +1 / L
    @inbounds k_g_l[1 , 7 ] = -1 / L
    @inbounds k_g_l[2 , 2 ] = +6 / (5 * L)
    @inbounds k_g_l[2 , 6 ] = +1 / 10
    @inbounds k_g_l[2 , 8 ] = -6 / (5 * L)
    @inbounds k_g_l[2 , 12] = +1 / 10
    @inbounds k_g_l[3 , 3 ] = +6 / (5 * L)
    @inbounds k_g_l[3 , 5 ] = -1 / 10
    @inbounds k_g_l[3 , 9 ] = -6 / (5 * L)
    @inbounds k_g_l[3 , 11] = -1 / 10
    @inbounds k_g_l[4 , 4 ] = +(I_zz + I_yy) / (A * L)
    @inbounds k_g_l[4 , 10] = -(I_zz + I_yy) / (A * L)
    @inbounds k_g_l[5 , 5 ] = +2 * L / 15
    @inbounds k_g_l[5 , 9 ] = +1 / 10
    @inbounds k_g_l[5 , 11] = -L / 30
    @inbounds k_g_l[6 , 6 ] = +2 * L / 15
    @inbounds k_g_l[6 , 8 ] = -1 / 10
    @inbounds k_g_l[6 , 12] = -L / 30
    @inbounds k_g_l[7 , 7 ] = +1 / L
    @inbounds k_g_l[8 , 8 ] = +6 / (5 * L)
    @inbounds k_g_l[8 , 12] = -1 / 10
    @inbounds k_g_l[9 , 9 ] = +6 / (5 * L)
    @inbounds k_g_l[9 , 11] = +1 / 10
    @inbounds k_g_l[10, 10] = +(I_zz + I_yy) / (A * L)
    @inbounds k_g_l[11, 11] = +2 * 1 * L / 15
    @inbounds k_g_l[12, 12] = +2 * 1 * L / 15

    # Compute the components of the element geometric stiffness matrix in its lower triangular part:
    for i in 1:12, j in (i + 1):12
        @inbounds k_g_l[j, i] = k_g_l[i, j]
    end

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, k_g_l, k_g_l)

    # Return the element geometric stiffness matrix:
    return k_g_l
end

function _compute_p_l(
    q_x::Real, q_y::Real, q_z::Real,
    L::Real)
    # Compute the element fixed-end force vector in the local coordinate system:
    p_l = [
        -q_x * L / 2       # F_x_i
        -q_y * L / 2       # F_y_i
        -q_z * L / 2       # F_z_i
        0                  # M_x_i
        -q_z * L ^ 2 / 12  # M_y_i
        -q_y * L ^ 2 / 12  # M_z_i
        +q_x * L / 2       # F_x_j
        -q_y * L / 2       # F_y_j
        -q_z * L / 2       # F_z_j
        0                  # M_x_j
        +q_z * L ^ 2 / 12  # M_y_j
        +q_y * L ^ 2 / 12] # M_z_j

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, p_l, p_l)

    # Return the element fixed-end forces:
    return p_l
end