"""
    struct Element

A type representing an element in the model of a structure of interest.

$(FIELDS)
"""
struct Element{CTI<:Real, CTJ<:Real, MPT<:Real, SPT<:Real}
    "Unique identifier"
    ID             ::Int
    "Unique identifier of the node ``i`` of the element"
    node_i_ID      ::Int
    "Unique identifier of the node ``j`` of the element"
    node_j_ID      ::Int
    "Unique identifier of the material of the element"
    material_ID    ::Int
    "Unique identifier of the section of the element"
    section_ID     ::Int
    "DOF releases at the node ``i``"
    releases_i     ::Vector{Bool}
    "DOF releases at the node ``j``"
    releases_j     ::Vector{Bool}
    "Angle that defines the orientation of the element's local coordinate system, ``\\omega``"
    ω              ::Real
    "``x``-coordinate of the node ``i``, ``x_{i}``"
    x_i            ::CTI
    "``y``-coordinate of the node ``i``, ``y_{i}``"
    y_i            ::CTI
    "``z``-coordinate of the node ``i``, ``z_{i}``"
    z_i            ::CTI
    "``x``-coordinate of the node ``j``, ``x_{j}``"
    x_j            ::CTJ
    "``y``-coordinate of the node ``j``, ``y_{j}``"
    y_j            ::CTJ
    "``z``-coordinate of the node ``j``, ``z_{j}``"
    z_j            ::CTJ
    "Young's modulus, ``E``"
    E              ::MPT
    "Poisson's ratio, ``\\nu``"
    ν              ::MPT
    "Density, ``\\rho``"
    ρ              ::MPT
    "Cross-sectional area, ``A``"
    A              ::SPT
    "Moment of inertia about the local ``z``-axis, ``I_{zz}``"
    I_zz           ::SPT
    "Moment of inertia about the local ``y``-axis, ``I_{yy}``"
    I_yy           ::SPT
    "Polar moment of inertia about the local ``x``-axis, ``J``"
    J              ::SPT
    "Length of the element, ``L``"
    L              ::Real
    "Local-to-global sub-transformation matrix, ``\\gamma``"
    γ              ::Matrix{<:Real}
    "Local-to-global transformation matrix, ``\\Gamma``"
    Γ              ::BlockDiagonal{<:Real}
    "Elastic stiffness matrix in the local coordinate system of the element, ``k_{e, l}``"
    k_e_l          ::Matrix{<:Real}
    "Elastic stiffness matrix in the global coordinate system, ``k_{e, g}``"
    k_e_g          ::Matrix{<:Real}
    "Geometric stiffness matrix in the local coordinate system of the element, ``k_{g, l}``"
    k_g_l          ::Matrix{<:Real}
    "Geometric stiffness matrix in the global coordinate system, ``k_{g, g}``"
    k_g_g          ::Matrix{<:Real}

    function Element(ID::Int, 
        node_i_ID::Int, node_j_ID::Int, material_ID::Int, section_ID::Int,
        releases_i::Vector{Bool}, releases_j::Vector{Bool}, ω::Real,
        x_i::CTI, y_i::CTI, z_i::CTI, 
        x_j::CTJ, y_j::CTJ, z_j::CTJ,
        E::MPT, ν::MPT, ρ::MPT, 
        A::SPT, I_zz::SPT, I_yy::SPT, J::SPT) where {CTI<:Real, CTJ<:Real, MPT<:Real, SPT<:Real}
        # Compute the length of the element:
        L = _compute_L(x_i, y_i, z_i, x_j, y_j, z_j)

        # Compute the transformation matrix:
        γ, Γ = _compute_Γ(x_i, y_i, z_i, x_j, y_j, z_j, L, ω)

        # Compute the element elastic stiffness matrix in its local coordinate system:
        k_e_l = _compute_k_e_l(E, ν, A, I_zz, I_yy, J, L)

        # Transform the element elastic stiffness matrix to the global coordinate system:
        k_e_g = Γ' * k_e_l * Γ

        # Remove small values if any:
        map!(x -> abs(x) < 1E-12 ? 0 : x, k_e_g, k_e_g)

        # Compute the element geometric stiffness matrix in its local coordinate system:
        k_g_l = _compute_k_g_l(A, I_zz, I_yy, L)

        # Transform the element geometric stiffness matrix to the global coordinate system:
        k_g_g = Γ' * k_g_l * Γ

        # Remove small values if any:
        map!(x -> abs(x) < 1E-12 ? 0 : x, k_g_g, k_g_g)

        # Return the element:
        return new{CTI, CTJ, MPT, SPT}(ID, 
            node_i_ID, node_j_ID, material_ID, section_ID, 
            releases_i, releases_j, ω, 
            x_i, y_i, z_i, 
            x_j, y_j, z_j, 
            E, ν, ρ, 
            A, I_zz, I_yy, J, 
            L, γ, Γ, k_e_l, k_e_g, k_g_l, k_g_g)
    end
end

function _compute_L(
    x_i::CTI, y_i::CTI, z_i::CTI,
    x_j::CTJ, y_j::CTJ, z_j::CTJ) where {CTI<:Real, CTJ<:Real}
    L = sqrt((x_j - x_i)^2 + (y_j - y_i)^2 + (z_j - z_i)^2)

    return L
end

function _compute_Γ(
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
        +c_χ * c_ρ                      +s_χ          -c_χ * s_ρ                  ;
        +s_ω * s_ρ - c_ω * s_χ * c_ρ    +c_ω * c_χ    +c_ω * s_χ * s_ρ + s_ω * c_ρ;
        +c_ω * s_ρ + s_ω * s_χ * c_ρ    -s_ω * c_χ    -s_ω * s_χ * s_ρ + c_ω * c_ρ]

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, γ, γ)

    # Preallocate the transformation matrix and fill it:
    Γ = BlockDiagonal([γ, γ, γ, γ])

    # Return the transformation matrices:
    return γ, Γ
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
    # map!(x -> abs(x) < 1E-12 ? 0 : x, k_g_l, k_g_l)

    # Return the element geometric stiffness matrix:
    return k_g_l
end

function _compute_p_l(
    q_x::Real, q_y::Real, q_z::Real,
    L::Real)
    # Compute the element fixed-end force vector in the local coordinate system:
    p_l = [
        -q_x * L / 2     ; # F_x_i
        -q_y * L / 2     ; # F_y_i
        -q_z * L / 2     ; # F_z_i
        0                ; # M_x_i
        -q_z * L ^ 2 / 12; # M_y_i
        -q_y * L ^ 2 / 12; # M_z_i
        +q_x * L / 2     ; # F_x_j
        -q_y * L / 2     ; # F_y_j
        -q_z * L / 2     ; # F_z_j
        0                ; # M_x_j
        +q_z * L ^ 2 / 12; # M_y_j
        +q_y * L ^ 2 / 12] # M_z_j

    # Remove small values if any:
    map!(x -> abs(x) < 1E-12 ? 0 : x, p_l, p_l)

    # Return the element fixed-end forces:
    return p_l
end