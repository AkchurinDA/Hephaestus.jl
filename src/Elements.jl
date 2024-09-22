"""
    struct Element

A type that represents an element in the finite element model of a structure.

This type should never be called directly by the user.

# Fields
$(FIELDS)
"""
struct Element{CTI<:Real, CTJ<:Real, MPT<:Real, SPT<:Real}
    # MAIN ELEMENT INFORMATION:
    "Unique identifier of the element provided by the user"
    ID          ::Int
    "Unique identifier of the node ``i`` of the element provided by the user"
    node_i_ID   ::Int
    "Unique identifier of the node ``j`` of the element provided by the user"
    node_j_ID   ::Int
    "Unique identifier of the material of the element provided by the user"
    material_ID ::Int
    "Unique identifier of the section of the element provided by the user"
    section_ID  ::Int
    "``x``-coordinate of the node ``i`` of the element, ``x_{i}``"
    x_i         ::CTI
    "``y``-coordinate of the node ``i`` of the element, ``y_{i}``"
    y_i         ::CTI
    "``z``-coordinate of the node ``i`` of the element, ``z_{i}``"
    z_i         ::CTI
    "``x``-coordinate of the node ``j`` of the element, ``x_{j}``"
    x_j         ::CTJ
    "``y``-coordinate of the node ``j`` of the element, ``y_{j}``"
    y_j         ::CTJ
    "``z``-coordinate of the node ``j`` of the element, ``z_{j}``"
    z_j         ::CTJ
    "Young's modulus of the material of the element, ``E``"
    E           ::MPT
    "Poisson's ratio of the material of the element, ``\\nu``"
    ν           ::MPT
    "Density of the material of the element, ``\\rho``"
    ρ           ::MPT
    "Cross-sectional area of the section of the element, ``A``"
    A           ::SPT
    "Moment of inertia about the local ``z``-axis of the section of the element, ``I_{zz}``"
    I_zz        ::SPT
    "Moment of inertia about the local ``y``-axis of the section of the element, ``I_{yy}``"
    I_yy        ::SPT
    "Polar moment of inertia of the section of the element, ``J``"
    J           ::SPT
    "Angle that defines the orientation of the local coordinate system of the element, ``\\omega``"
    ω           ::Real
    "DOF releases at the node ``i`` of the element"
    releases_i  ::Vector{Bool}
    "DOF releases at the node ``j`` of the element"
    releases_j  ::Vector{Bool}

    # ADDITIONAL ELEMENT INFORMATION THAT CAN BE PRECOMPUTED:
    "Length of the element, ``L``"
    L           ::Real
    "Global-to-local sub-transformation matrix of the element, ``[\\gamma]``"
    γ           ::Matrix{<:Real}
    "Global-to-local transformation matrix of the element, ``[T]``"
    T           ::Matrix{<:Real}
    "Element's elastic stiffness matrix in its local coordinate system, ``[k_{e, l}]``"
    k_e_l       ::Matrix{<:Real}
    "Element's geometric stiffness matrix in its local coordinate system, ``[k_{g, l}]``"
    k_g_l       ::Matrix{<:Real}
    "Element's elastic stiffness matrix in its global coordinate system, ``[k_{e, g}]``"
    k_e_g       ::Matrix{<:Real}
    "Element's geometric stiffness matrix in its global coordinate system, ``[k_{g, g}]``"
    k_g_g       ::Matrix{<:Real}
    # "Condensed element's elastic stiffness matrix in its local coordinate system, ``[k_{e, l, c}]``"
    # k_e_l_c     ::Matrix{<:Real}
    # "Condensed element's geometric stiffness matrix in its local coordinate system, ``[k_{g, l, c}]``"
    # k_g_l_c     ::Matrix{<:Real}
end

"""
    _compute_T()

This function computes the global-to-local transformation matrix of an element.
"""
function _compute_T(
    x_i::CTI, y_i::CTI, z_i::CTI,
    x_j::CTJ, y_j::CTJ, z_j::CTJ,
    L::Real,
    ω::Real) where {CTI<:Real, CTJ<:Real}
    # Compute the rotation angles:
        # ρ - 1st rotation about local z-axis
        # χ - 2nd rotation about local y-axis
        # ω - 3rd rotation about local x-axis
    ρ = -atan(z_j - z_i, x_j - x_i)
    χ = π / 2 - acos((y_j - y_i) / L)

    # Construct the sub-transformation matrix:
        # γ = γ_ω * γ_χ * γ_ρ
    s_ρ, c_ρ = sincos(ρ)
    s_χ, c_χ = sincos(χ)
    s_ω, c_ω = sincos(ω)
    γ = [
        +c_χ * c_ρ                      +s_χ          -c_χ * s_ρ                  ;
        +s_ω * s_ρ - c_ω * s_χ * c_ρ    +c_ω * c_χ    +c_ω * s_χ * s_ρ + s_ω * c_ρ;
        +c_ω * s_ρ + s_ω * s_χ * c_ρ    -s_ω * c_χ    -s_ω * s_χ * s_ρ + c_ω * c_ρ]

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, γ, γ)

    # Preallocate the transformation matrix and fill it:
    T = zeros(eltype(γ), 12, 12)
    @inbounds T[1:3  , 1:3  ] = γ
    @inbounds T[4:6  , 4:6  ] = γ
    @inbounds T[7:9  , 7:9  ] = γ
    @inbounds T[10:12, 10:12] = γ

    # Return the transformation matrix:
    return γ, T
end

"""
    _compute_k_e_l()

This function computes the element elastic stiffness matrix in its local coordinate system.
"""
function _compute_k_e_l(
    E::MPT, ν::MPT,
    A::SPT, I_zz::SPT, I_yy::SPT, J::SPT,
    L::ELT) where {MPT<:Real, SPT<:Real, ELT<:Real}
    # Preallocate:
    T     = float(promote_type(MPT, SPT, ELT))
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
    
    # Return the element stiffness matrix:
    return k_e_l
end

"""
    _compute_k_g_l()

This function computes the element geometric stiffness matrix in its local coordinate system (without the axial force term ``P``).
"""
function _compute_k_g_l(
    A::SPT, I_zz::SPT, I_yy::SPT,
    L::ELT) where {SPT<:Real, ELT<:Real}
    # # Preallocate:
    T     = float(promote_type(SPT, ELT))
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

"""
    _compute_p_l()

This function computes the element fixed-end forces due to applied distributed loads.
"""
function _compute_p_l(
    w_x::DLTX, w_y::DLTY, w_z::DLTZ, 
    L::ELT) where {DLTX<:Real, DLTY<:Real, DLTZ<:Real, ELT<:Real}
    # Compute the element fixed-end forces:
    p_l = [
        -w_x * L / 2     ; # F_x_i
        -w_y * L / 2     ; # F_y_i
        -w_z * L / 2     ; # F_z_i
        0                ; # M_x_i
        -w_z * L ^ 2 / 12; # M_y_i
        -w_y * L ^ 2 / 12; # M_z_i
        +w_x * L / 2     ; # F_x_j
        -w_y * L / 2     ; # F_y_j
        -w_z * L / 2     ; # F_z_j
        0                ; # M_x_j
        +w_z * L ^ 2 / 12; # M_y_j
        +w_y * L ^ 2 / 12] # M_z_j

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, p_l, p_l)

    # Return the element fixed-end forces:
    return p_l
end