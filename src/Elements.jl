abstract type AbstractElement end

"""
    struct Element

A type representing an element in the finite element model.

$(FIELDS)
"""
struct Element{
    NCTI  <:Real, # NCTI   - (N)ode (C)oordinate (T)ype (I)
    NCTJ  <:Real, # NCTJ   - (N)ode (C)oordinate (T)ype (J)
    MPT   <:Real, # MPT    - (M)aterial (P)roperty (T)ype
    SPT   <:Real, # SPT    - (S)ection (P)roperty (T)ype
    LTGTMT<:Real, # LTGTMT - (L)ocal (T)o (G)lobal (T)ransformation (M)atrix (T)ype
    EGESMT<:Real, # EGESMT - (E)lement (G)lobal (E)lastic (S)tiffness (M)atrix (T)ype
    EGGSMT<:Real, # EGGSMT - (E)lement (G)lobal (G)eometric (S)tiffness (M)atrix (T)ype
    EGMMT <:Real  # EGMMT  - (E)lement (G)lobal (M)ass (M)atrix (T)ype
    } <: AbstractElement
    "Unique identifier of the node (``i``)"
    node_i_tag  ::Int
    "``x``-coordinate of the node (``i``)"
    x_i         ::NCTI
    "``y``-coordinate of the node (``i``)"
    y_i         ::NCTI
    "``z``-coordinate of the node (``i``)"
    z_i         ::NCTI
    "Unique identifier of the node (``j``)"
    node_j_tag  ::Int
    "``x``-coordinate of the node (``j``)"
    x_j         ::NCTJ
    "``y``-coordinate of the node (``j``)"
    y_j         ::NCTJ
    "``z``-coordinate of the node (``j``)"
    z_j         ::NCTJ
    "Unique identifier of the material"
    material_tag::Int
    "Young's modulus"
    E           ::MPT
    "Poisson's ratio"
    ν           ::MPT
    "Density"
    ρ           ::MPT
    "Unique identifier of the section"
    section_tag ::Int
    "Cross-sectional area"
    A           ::SPT
    "Moment of inertia about the local ``z``-axis"
    I_zz        ::SPT
    "Moment of inertia about the local ``x``-axis"
    I_yy        ::SPT
    "Polar moment of inertia"
    J           ::SPT
    "Angle that defines the orientation of the local coordinate system of the element"
    ω           ::Real
    "DOF releases at the node (``i``)"
    releases_i  ::Vector{Bool}
    "DOF releases at the node (``j``)"
    releases_j  ::Vector{Bool}
    "Length of the element"
    L           ::Real
    "Local-to-global sub-transformation matrix"
    γ           ::AbstractMatrix{LTGTMT}
    "Local-to-global transformation matrix"
    Γ           ::AbstractMatrix{LTGTMT}
    "Element elastic stiffness matrix in the local coordinate system"
    k_e_l       ::AbstractMatrix{<:Real}
    "Element elastic stiffness matrix in the global coordinate system"
    k_e_g       ::AbstractMatrix{EGESMT}
    "Element geometric stiffness matrix in the local coordinate system"
    k_g_l       ::AbstractMatrix{<:Real}
    "Element geometric stiffness matrix in the global coordinate system"
    k_g_g       ::AbstractMatrix{EGGSMT}
    "Element mass matrix in the local coordinate system"
    m_l         ::AbstractMatrix{<:Real}
    "Element mass matrix in the global coordinate system"
    m_g         ::AbstractMatrix{EGMMT }

    function Element(
        node_i_tag::Int, x_i::NCTI, y_i::NCTI, z_i::NCTI,
        node_j_tag::Int, x_j::NCTJ, y_j::NCTJ, z_j::NCTJ,
        material_tag::Int, E::MPT, ν::MPT, ρ::MPT,
        section_tag::Int, A::SPT, I_zz::SPT, I_yy::SPT, J::SPT,
        ω::Real,
        releases_i::Vector{Bool}, releases_j::Vector{Bool}) where {NCTI<:Real, NCTJ<:Real, MPT<:Real, SPT<:Real}
        L = _compute_L(x_i, y_i, z_i, x_j, y_j, z_j)

        γ = _compute_γ(x_i, y_i, z_i, x_j, y_j, z_j, ω, L)
        Γ = _compute_Γ(γ)

        k_e_l = _compute_k_e_l(E, ν, A, I_zz, I_yy, J, L)
        k_e_g = transpose(Γ) * k_e_l * Γ

        k_g_l = _compute_k_g_l(A, J, L)
        k_g_g = transpose(Γ) * k_g_l * Γ

        m_l = _compute_m_l(ρ, A, J, L)
        m_g = transpose(Γ) * m_l * Γ

        LTGTMT = eltype(Γ)
        EGESMT = eltype(k_e_g)
        EGGSMT = eltype(k_g_g)
        EGMMT  = eltype(m_g)

        return new{NCTI, NCTJ, MPT, SPT, LTGTMT, EGESMT, EGGSMT, EGMMT}(
            node_i_tag, x_i, y_i, z_i,
            node_j_tag, x_j, y_j, z_j,
            material_tag, E, ν, ρ,
            section_tag, A, I_zz, I_yy, J,
            ω,
            releases_i, releases_j,
            L,
            γ, Γ,
            k_e_l, k_e_g, k_g_l, k_g_g, m_l, m_g)
    end
end

function _compute_L(
    x_i::NCTI, y_i::NCTI, z_i::NCTI,
    x_j::NCTJ, y_j::NCTJ, z_j::NCTJ) where {NCTI<:Real, NCTJ<:Real}
    L = sqrt((x_j - x_i)^2 + (y_j - y_i)^2 + (z_j - z_i)^2)

    return L
end

function _compute_γ(
    x_i::NCTI, y_i::NCTI, z_i::NCTI,
    x_j::NCTJ, y_j::NCTJ, z_j::NCTJ,
    ω::Real,
    L::EPT) where {NCTI<:Real, NCTJ<:Real, EPT<:Real}
    Δx = x_j - x_i
    Δy = y_j - y_i
    Δz = z_j - z_i

    ρ = -atan(Δz / Δx)
    χ = π / 2 - acos(Δy / L)

    s_ρ, c_ρ = sincos(ρ)
    s_χ, c_χ = sincos(χ)
    s_ω, c_ω = sincos(ω)

    γ = [
        +c_χ * c_ρ                      +s_χ          -c_χ * s_ρ                  ;
        -c_ω * s_χ * c_ρ + s_ω * s_ρ    +c_ω * c_χ    +c_ω * s_χ * s_ρ + s_ω * c_ρ;
        +s_ω * s_χ * c_ρ + c_ω * s_ρ    -s_ω * c_χ    -s_ω * s_χ * s_ρ + c_ω * c_ρ]

    return γ
end

function _compute_Γ(γ::AbstractMatrix{T}) where {T<:Real}
    Γ = zeros(T, 12, 12)

    Γ[1:3  , 1:3  ] .= γ
    Γ[4:6  , 4:6  ] .= γ
    Γ[7:9  , 7:9  ] .= γ
    Γ[10:12, 10:12] .= γ

    return Γ
end

function _compute_k_e_l(
    E::MPT, ν::MPT, 
    A::SPT, I_zz::SPT, I_yy::SPT, J::SPT, 
    L::EPT) where {MPT<:Real, SPT<:Real, EPT<:Real}
    T     = float(promote_type(MPT, SPT, EPT))
    k_e_l = zeros(T, 12, 12)

    __compute_k_e_l!(k_e_l, E, ν, A, I_zz, I_yy, J, L)

    return k_e_l
end

function __compute_k_e_l!(k_e_l::AbstractMatrix{ELESMT},
    E::MPT, ν::MPT, 
    A::SPT, I_zz::SPT, I_yy::SPT, J::SPT, 
    L::EPT) where {ELESMT<:Real, MPT<:Real, SPT<:Real, EPT<:Real}
    G = E / (2 * (1 + ν))

    @inbounds k_e_l[1 , 1 ] = +E * A / L
    @inbounds k_e_l[1 , 7 ] = -E * A / L
    @inbounds k_e_l[2 , 2 ] = +12 * E * I_zz / L ^ 3
    @inbounds k_e_l[2 , 6 ] = +6  * E * I_zz / L ^ 2
    @inbounds k_e_l[2 , 8 ] = -12 * E * I_zz / L ^ 3
    @inbounds k_e_l[2 , 12] = +6  * E * I_zz / L ^ 2
    @inbounds k_e_l[3 , 3 ] = +12 * E * I_yy / L ^ 3
    @inbounds k_e_l[3 , 5 ] = -6  * E * I_yy / L ^ 2
    @inbounds k_e_l[3 , 9 ] = -12 * E * I_yy / L ^ 3
    @inbounds k_e_l[3 , 11] = -6  * E * I_yy / L ^ 2
    @inbounds k_e_l[4 , 4 ] = +G * J / L
    @inbounds k_e_l[4 , 10] = -G * J / L
    @inbounds k_e_l[5 , 5 ] = +4  * E * I_yy / L
    @inbounds k_e_l[5 , 9 ] = +6  * E * I_yy / L ^ 2
    @inbounds k_e_l[5 , 11] = +2  * E * I_yy / L
    @inbounds k_e_l[6 , 6 ] = +4  * E * I_zz / L
    @inbounds k_e_l[6 , 8 ] = -6  * E * I_zz / L ^ 2
    @inbounds k_e_l[6 , 12] = +2  * E * I_zz / L
    @inbounds k_e_l[7 , 7 ] = +E * A / L
    @inbounds k_e_l[8 , 8 ] = +12 * E * I_zz / L ^ 3
    @inbounds k_e_l[8 , 12] = -6  * E * I_zz / L ^ 2
    @inbounds k_e_l[9 , 9 ] = +12 * E * I_yy / L ^ 3
    @inbounds k_e_l[9 , 11] = +6  * E * I_yy / L ^ 2
    @inbounds k_e_l[10, 10] = +G * J / L
    @inbounds k_e_l[11, 11] = +4  * E * I_yy / L
    @inbounds k_e_l[12, 12] = +4  * E * I_zz / L

    for i in 1:12, j in (i + 1):12
        @inbounds k_e_l[j, i] = k_e_l[i, j]
    end

    return k_e_l
end

function _compute_k_g_l(
    A::SPT, J::SPT,
    L::EPT) where {SPT<:Real, EPT<:Real}
    T     = float(promote_type(SPT, EPT))
    k_g_l = zeros(T, 12, 12)

    __compute_k_g_l!(k_g_l, A, J, L)

    return k_g_l
end

function __compute_k_g_l!(k_g_l::AbstractMatrix{ELGSMT},
    A::SPT, J::SPT,
    L::EPT) where {ELGSMT<:Real, SPT<:Real, EPT<:Real}
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
    @inbounds k_g_l[4 , 4 ] = +J / (A * L)
    @inbounds k_g_l[4 , 10] = -J / (A * L)
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
    @inbounds k_g_l[10, 10] = +J / (A * L)
    @inbounds k_g_l[11, 11] = +2 * L / 15
    @inbounds k_g_l[12, 12] = +2 * L / 15

    # Compute the components of the element geometric stiffness matrix in its lower triangular part:
    for i in 1:12, j in (i + 1):12
        @inbounds k_g_l[j, i] = k_g_l[i, j]
    end

    return k_g_l
end


function _compute_m_l(
    ρ::MPT,
    A::SPT, J::SPT,
    L::EPT) where {MPT<:Real, SPT<:Real, EPT<:Real}
    T   = float(promote_type(MPT, SPT, EPT))
    m_l = zeros(T, 12, 12)

    __compute_m_l!(m_l, ρ, A, J, L)

    return m_l
end

function __compute_m_l!(m_l::AbstractMatrix{<:Real},
    ρ::MPT,
    A::SPT, J::SPT,
    L::EPT) where {MPT<:Real, SPT<:Real, EPT<:Real}
    @inbounds m_l[1 , 1 ] = +140
    @inbounds m_l[1 , 7 ] = +70
    @inbounds m_l[2 , 2 ] = +156
    @inbounds m_l[2 , 6 ] = +22 * L
    @inbounds m_l[2 , 8 ] = +54
    @inbounds m_l[2 , 12] = -13 * L
    @inbounds m_l[3 , 3 ] = +156
    @inbounds m_l[3 , 5 ] = -22 * L
    @inbounds m_l[3 , 9 ] = +54 
    @inbounds m_l[3 , 11] = +13 * L
    @inbounds m_l[4 , 4 ] = +140 * J / A
    @inbounds m_l[4 , 10] = +70  * J / A
    @inbounds m_l[5 , 5 ] = +4 * L ^ 2 
    @inbounds m_l[5 , 9 ] = -13 * L
    @inbounds m_l[5 , 11] = -3 * L ^ 2
    @inbounds m_l[6 , 6 ] = +4 * L ^ 2
    @inbounds m_l[6 , 8 ] = +13 * L
    @inbounds m_l[6 , 12] = -3 * L ^ 2
    @inbounds m_l[7 , 7 ] = +140
    @inbounds m_l[8 , 8 ] = +156
    @inbounds m_l[8 , 12] = -22 * L
    @inbounds m_l[9 , 9 ] = +156
    @inbounds m_l[9 , 11] = +22 * L
    @inbounds m_l[10, 10] = +140 * J / A
    @inbounds m_l[11, 11] = +4 * L ^ 2
    @inbounds m_l[12, 12] = +4 * L ^ 2

    for i in 1:12, j in (i + 1):12
        @inbounds m_l[j, i] = m_l[i, j]
    end

    m_l .*= (ρ * A * L) / 420

    return m_l
end


function _compute_p_l(
    q_x::Real, q_y::Real, q_z::Real,
    L::Real)
    # TODO: Add distributed moments

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

    # Return the element fixed-end forces:
    return p_l
end

_get_k_e_g_type(::Element{NCTI, NCTJ, MPT, SPT, LTGTMT, EGESMT, EGGSMT, EGMMT}) where {
    NCTI, NCTJ, 
    MPT, 
    SPT, 
    LTGTMT, EGESMT, EGGSMT, EGMMT} = EGESMT
_get_k_g_g_type(::Element{NCTI, NCTJ, MPT, SPT, LTGTMT, EGESMT, EGGSMT, EGMMT}) where {
    NCTI, NCTJ, 
    MPT, 
    SPT, 
    LTGTMT, EGESMT, EGGSMT, EGMMT} = EGGSMT
_get_m_g_type(::Element{NCTI, NCTJ, MPT, SPT, LTGTMT, EGESMT, EGGSMT, EGMMT}) where {
    NCTI, NCTJ, 
    MPT, 
    SPT, 
    LTGTMT, EGESMT, EGGSMT, EGMMT} = EGMMT