"""
    struct Element

A type representing an element in the finite element model.

$(FIELDS)
"""
struct Element{
    NIT <: Real, # (N)ode (i) (T)ype
    NJT <: Real, # (N)ode (j) (T)ype
    ST  <: Real, # (S)ection (T)ype
    MT  <: Real, # (M)aterial (T)ype
    OAT <: Real, # (O)rientation (A)ngle (T)ype
    ET  <: Real} # (E)lement (T)ype
    "Identification tag"
    ID        ::Int
    "Node (``i``) of the element"
    node_i    ::Node{NIT}
    "Node (``j``) of the element"
    node_j    ::Node{NJT}
    "Section attached to the element"
    section   ::Section{ST}
    "Material attached to the element"
    material  ::Material{MT}
    "Orientation angle of the element"
    ω         ::OAT
    "End moment releases at node (``i``) of the element"
    releases_i::Vector{Bool}
    "End moment releases at node (``j``) of the element"
    releases_j::Vector{Bool}

    function Element(ID, 
        node_i    ::Node{NIT},
        node_j    ::Node{NJT},
        section   ::Section{ST},
        material  ::Material{MT},
        ω         ::OAT,
        releases_i::Vector{<:Bool},
        releases_j::Vector{<:Bool})::Element where {
        NIT <: Real,
        NJT <: Real,
        ST  <: Real,
        MT  <: Real,
        OAT <: Real}
        # Promote the types:
        ET = float(promote_type(NIT, NJT, ST, MT, OAT))

        # Return the element:
        return new{NIT, NJT, ST, MT, OAT, ET}(ID, node_i, node_j, section, material, ω, releases_i, releases_j)
    end
end

gettype(::Element{NIT, NJT, ST, MT, OAT, ET}) where {NIT, NJT, ST, MT, OAT, ET} = ET

mutable struct ElementState
    "Identification tag"
    ID   ::Int
    "Current length of the element"
    L    ::Real
    "Current orientation angle of the element at node (``i``)"
    ω_i  ::Real
    "Current orientation angle of the element at node (``j``)"
    ω_j  ::Real
    "Current element global-to-local transformation matrix"
    Γ    ::AbstractMatrix{<:Real}
    "Current element internal force vector in its local coordinate system"
    q    ::AbstractVector{<:Real}
    "Current element elastic stiffness matrix in its local coordinate system"
    k_e_l::AbstractMatrix{<:Real}
    "Current element elastic stiffness matrix in the global coordinate system"
    k_e_g::AbstractMatrix{<:Real}
    "Current element geometric stiffness matrix in its local coordinate system"
    k_g_l::AbstractMatrix{<:Real}
    "Current element geometric stiffness matrix in the global coordinate system"
    k_g_g::AbstractMatrix{<:Real}
    modified::Bool
end

function compute_γ(ρ::Real, χ::Real, ω::Real)::AbstractMatrix{<:Real}
    # Compute the rotation matrix about the y-axis:
    s_ρ, c_ρ = sincos(ρ)
    γ_ρ = [+c_ρ 0 -s_ρ; 0 +1 0; +s_ρ 0 +c_ρ]

    # Compute the rotation matrix about the z-axis:
    s_χ, c_χ = sincos(χ)
    γ_χ = [+c_χ +s_χ 0; -s_χ +c_χ 0; 0 0 +1]

    # Compute the rotation matrix about the x-axis:
    s_ω, c_ω = sincos(ω)
    γ_ω = [+1 0 0; 0 +c_ω +s_ω; 0 -s_ω +c_ω]

    # Compute the element global-to-local subtransformation matrix:
    γ = γ_ω * γ_χ * γ_ρ

    # Return the element global-to-local subtransformation matrix:
    return γ
end

function compute_Γ(
    x_i::Real, y_i::Real, z_i::Real, ω_i::Real,
    x_j::Real, y_j::Real, z_j::Real, ω_j::Real,
    L::Real)::AbstractMatrix{<:Real}
    # Compute the element length projections:
    Δx = x_j - x_i
    Δy = y_j - y_i
    Δz = z_j - z_i

    # Compute the element orientation angles:
    ρ = -atan(Δz, Δx)
    χ = π / 2 - acos(Δy / L)

    # Compute the element global-to-local subtransformation matrix:
    γ_i = compute_γ(ρ, χ, ω_i) # At node (i)
    γ_j = compute_γ(ρ, χ, ω_j) # At node (j)

    # Compute the element global-to-local transformation matrix:
    T = promote_type(eltype(γ_i), eltype(γ_j))
    Γ = zeros(T, 12, 12)
    Γ[  1:3,   1:3] .= γ_i
    Γ[  4:6,   4:6] .= γ_i
    Γ[  7:9,   7:9] .= γ_j
    Γ[10:12, 10:12] .= γ_j

    # Return the element global-to-local transformation matrix:
    return Γ
end

function compute_k_e_l!(k_e_l::AbstractMatrix{<:Real}, element::Element, L::Real)::AbstractMatrix{<:Real}
    # Extract the material properties:
    material = element.material
    E, ν = material.E, material.ν

    # Extract the section properties:
    section = element.section
    A, I_zz, I_yy, J = section.A, section.I_zz, section.I_yy, section.J

    # Compute the element elastic stiffness matrix:
    compute_k_e_l!(k_e_l, E, ν, A, I_zz, I_yy, J, L)

    # Return the element elastic stiffness matrix:
    return k_e_l
end

function compute_k_e_l!(k_e_l::AbstractMatrix{<:Real},
    E::Real, ν::Real, 
    A::Real, I_zz::Real, I_yy::Real, J::Real, 
    L::Real)::AbstractMatrix{<:Real}
    # Compute the shear modulus:
    G = E / (2 * (1 + ν))

    # Construct the upper triangular part of the element elastic stiffness matrix:
    @inbounds k_e_l[ 1,  1] = +E * A / L
    @inbounds k_e_l[ 1,  7] = -E * A / L
    @inbounds k_e_l[ 2,  2] = +12 * E * I_zz / L ^ 3
    @inbounds k_e_l[ 2,  6] = +6  * E * I_zz / L ^ 2
    @inbounds k_e_l[ 2,  8] = -12 * E * I_zz / L ^ 3
    @inbounds k_e_l[ 2, 12] = +6  * E * I_zz / L ^ 2
    @inbounds k_e_l[ 3,  3] = +12 * E * I_yy / L ^ 3
    @inbounds k_e_l[ 3,  5] = -6  * E * I_yy / L ^ 2
    @inbounds k_e_l[ 3,  9] = -12 * E * I_yy / L ^ 3
    @inbounds k_e_l[ 3, 11] = -6  * E * I_yy / L ^ 2
    @inbounds k_e_l[ 4,  4] = +G * J / L
    @inbounds k_e_l[ 4, 10] = -G * J / L
    @inbounds k_e_l[ 5,  5] = +4  * E * I_yy / L
    @inbounds k_e_l[ 5,  9] = +6  * E * I_yy / L ^ 2
    @inbounds k_e_l[ 5, 11] = +2  * E * I_yy / L
    @inbounds k_e_l[ 6,  6] = +4  * E * I_zz / L
    @inbounds k_e_l[ 6,  8] = -6  * E * I_zz / L ^ 2
    @inbounds k_e_l[ 6, 12] = +2  * E * I_zz / L
    @inbounds k_e_l[ 7,  7] = +E * A / L
    @inbounds k_e_l[ 8,  8] = +12 * E * I_zz / L ^ 3
    @inbounds k_e_l[ 8, 12] = -6  * E * I_zz / L ^ 2
    @inbounds k_e_l[ 9,  9] = +12 * E * I_yy / L ^ 3
    @inbounds k_e_l[ 9, 11] = +6  * E * I_yy / L ^ 2
    @inbounds k_e_l[10, 10] = +G * J / L
    @inbounds k_e_l[11, 11] = +4  * E * I_yy / L
    @inbounds k_e_l[12, 12] = +4  * E * I_zz / L

    # Construct the lower triangular part of the element elastic stiffness matrix:
    for i in 1:12, j in (i + 1):12
        @inbounds k_e_l[j, i] = k_e_l[i, j]
    end

    # Return the element elastic stiffness matrix:
    return k_e_l
end

function compute_k_g_l!(k_g_l::AbstractMatrix{<:Real}, element::Element, L::Real, N::Real)::AbstractMatrix{<:Real}
    # Extract the section properties:
    section = element.section
    A, J = section.A, section.J

    # Compute the element elastic stiffness matrix:
    compute_k_g_l!(k_g_l, A, J, L, N)

    # Return the element elastic stiffness matrix:
    return k_g_l
end

function compute_k_g_l!(k_g_l::AbstractMatrix{<:Real},
    A::Real, J::Real,
    L::Real,
    N::Real)::AbstractMatrix{<:Real}
    # Compute the upper triangular part of the element geometric stiffness matrix:
    @inbounds k_g_l[ 1,  1] = +1 / L
    @inbounds k_g_l[ 1,  7] = -1 / L
    @inbounds k_g_l[ 2,  2] = +6 / (5 * L)
    @inbounds k_g_l[ 2,  6] = +1 / 10
    @inbounds k_g_l[ 2,  8] = -6 / (5 * L)
    @inbounds k_g_l[ 2, 12] = +1 / 10
    @inbounds k_g_l[ 3,  3] = +6 / (5 * L)
    @inbounds k_g_l[ 3,  5] = -1 / 10
    @inbounds k_g_l[ 3,  9] = -6 / (5 * L)
    @inbounds k_g_l[ 3, 11] = -1 / 10
    @inbounds k_g_l[ 4,  4] = +J / (A * L)
    @inbounds k_g_l[ 4, 10] = -J / (A * L)
    @inbounds k_g_l[ 5,  5] = +2 * L / 15
    @inbounds k_g_l[ 5,  9] = +1 / 10
    @inbounds k_g_l[ 5, 11] = -L / 30
    @inbounds k_g_l[ 6,  6] = +2 * L / 15
    @inbounds k_g_l[ 6,  8] = -1 / 10
    @inbounds k_g_l[ 6, 12] = -L / 30
    @inbounds k_g_l[ 7,  7] = +1 / L
    @inbounds k_g_l[ 8,  8] = +6 / (5 * L)
    @inbounds k_g_l[ 8, 12] = -1 / 10
    @inbounds k_g_l[ 9,  9] = +6 / (5 * L)
    @inbounds k_g_l[ 9, 11] = +1 / 10
    @inbounds k_g_l[10, 10] = +J / (A * L)
    @inbounds k_g_l[11, 11] = +2 * L / 15
    @inbounds k_g_l[12, 12] = +2 * L / 15

    k_g_l .*= +N

    # Compute the lower triangular part of the element geometric stiffness matrix:
    for i in 1:12, j in (i + 1):12
        @inbounds k_g_l[j, i] = k_g_l[i, j]
    end

    # Return the element geometric stiffness matrix:
    return k_g_l
end

function compute_m_l!(m_l::AbstractMatrix{<:Real},
    ρ::Real, 
    A::Real, J::Real, 
    L::Real)::AbstractMatrix{<:Real}
    # Compute the upper triangular part of the element mass matrix:
    @inbounds m_l[ 1,  1] = +140
    @inbounds m_l[ 1,  7] = +70
    @inbounds m_l[ 2,  2] = +156
    @inbounds m_l[ 2,  6] = +22 * L
    @inbounds m_l[ 2,  8] = +54
    @inbounds m_l[ 2, 12] = -13 * L
    @inbounds m_l[ 3,  3] = +156
    @inbounds m_l[ 3,  5] = -22 * L
    @inbounds m_l[ 3,  9] = +54 
    @inbounds m_l[ 3, 11] = +13 * L
    @inbounds m_l[ 4,  4] = +140 * J / A
    @inbounds m_l[ 4, 10] = +70  * J / A
    @inbounds m_l[ 5,  5] = +4 * L ^ 2 
    @inbounds m_l[ 5,  9] = -13 * L
    @inbounds m_l[ 5, 11] = -3 * L ^ 2
    @inbounds m_l[ 6,  6] = +4 * L ^ 2
    @inbounds m_l[ 6,  8] = +13 * L
    @inbounds m_l[ 6, 12] = -3 * L ^ 2
    @inbounds m_l[ 7,  7] = +140
    @inbounds m_l[ 8,  8] = +156
    @inbounds m_l[ 8, 12] = -22 * L
    @inbounds m_l[ 9,  9] = +156
    @inbounds m_l[ 9, 11] = +22 * L
    @inbounds m_l[10, 10] = +140 * J / A
    @inbounds m_l[11, 11] = +4 * L ^ 2
    @inbounds m_l[12, 12] = +4 * L ^ 2

    m_l .*= +ρ * A * L / 420

    # Compute the lower triangular part of the element mass matrix:
    for i in 1:12, j in (i + 1):12
        @inbounds m_l[j, i] = m_l[i, j]
    end

    # Return the element mass matrix:
    return m_l
end

function condense!(m::AbstractMatrix{<:Real}, releases_i::Vector{Bool}, releases_j::Vector{Bool})::AbstractMatrix{<:Real}
    # Condense the matrix if end releases are present:
    for i in 1:3
        # Node (i):
        if releases_i[i]
            m[3 + i, :] .= 0
            m[:, 3 + i] .= 0
        end

        # Node (j):
        if releases_j[i]
            m[9 + i, :] .= 0
            m[:, 9 + i] .= 0
        end
    end

    # Return the condensed matrix:
    return m
end

function transform(m::AbstractMatrix{<:Real}, Γ::AbstractMatrix{<:Real})::AbstractMatrix{<:Real}
    # Transform the matrix to the global coordinate system:
    M = transpose(Γ) * m * Γ

    # Return the transformed matrix:
    return M
end

function initelementstate(element::Element)::ElementState
    # Compute the length of the element:
    L = sqrt(
        (element.node_j.x - element.node_i.x) ^ 2 + 
        (element.node_j.y - element.node_i.y) ^ 2 + 
        (element.node_j.z - element.node_i.z) ^ 2)

    # Compute the orientation angle of the element:
    ω_i = element.ω
    ω_j = element.ω

    # Compute the element global-to-local transformation matrix:
    Γ = compute_Γ(
        element.node_i.x, element.node_i.y, element.node_i.z, ω_i, 
        element.node_j.x, element.node_j.y, element.node_j.z, ω_j, 
        L)

    # Infer the type of the element stiffness matrices:
    T = gettype(element)

    # Initialize the element internal force vector:
    q = zeros(T, 12)

    # Initialize the element elastic stiffness matrix:
    k_e_l = zeros(T, 12, 12)
    compute_k_e_l!(k_e_l, element, L)
    condense!(k_e_l, element.releases_i, element.releases_j)
    k_e_g = transform(k_e_l, Γ)

    # Initialize the element geometric stiffness matrix:
    k_g_l = zeros(T, 12, 12)
    compute_k_g_l!(k_g_l, element, L, 0)
    condense!(k_g_l, element.releases_i, element.releases_j)
    k_g_g = transform(k_g_l, Γ)

    # Initialize the state of the element:
    elementstate = ElementState(element.ID, L, ω_i, ω_j, Γ, q, k_e_l, k_e_g, k_g_l, k_g_g, false)

    # Return the state of the element:
    return elementstate
end