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
    "Has the state been modified?"
    modified::Bool

    # Constructor:
    ElementState() = new()
end

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
    "Current state"
    state     ::ElementState

    # Constructor:
    function Element(ID,
        node_i    ::Node{NIT},
        node_j    ::Node{NJT},
        section   ::Section{ST},
        material  ::Material{MT},
        ω         ::OAT,
        releases_i::Vector{<:Bool},
        releases_j::Vector{<:Bool},
        state     ::ElementState)::Element where {
        NIT <: Real,
        NJT <: Real,
        ST  <: Real,
        MT  <: Real,
        OAT <: Real}
        # Promote the types:
        ET = float(promote_type(NIT, NJT, ST, MT, OAT))

        # Return the element:
        return new{NIT, NJT, ST, MT, OAT, ET}(ID, node_i, node_j, section, material, ω, releases_i, releases_j, state)
    end
end

function initstate!(element::Element)::Element
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
    element.state.ID       = element.ID
    element.state.L        = L
    element.state.ω_i      = ω_i
    element.state.ω_j      = ω_j
    element.state.Γ        = Γ
    element.state.q        = q
    element.state.k_e_l    = k_e_l
    element.state.k_e_g    = k_e_g
    element.state.k_g_l    = k_g_l
    element.state.k_g_g    = k_g_g
    element.state.modified = false

    # Return the state of the element:
    return element
end

function updatestate!(element::Element, δu_i::AbstractVector{T}, δu_j::AbstractVector{T})::Element where {T<:Real}
    # Extract the original nodal coordinates:
    x_i = element.node_i.x
    y_i = element.node_i.y
    z_i = element.node_i.z
    x_j = element.node_j.x
    y_j = element.node_j.y
    z_j = element.node_j.z

    # Extract the nodal coordinates:
    u_x_i = element.node_i.state.u_x
    u_y_i = element.node_i.state.u_y
    u_z_i = element.node_i.state.u_z
    u_x_j = element.node_j.state.u_x
    u_y_j = element.node_j.state.u_y
    u_z_j = element.node_j.state.u_z

    # Update the element length:
    element.state.L = sqrt(
        (x_j + u_x_i - x_i - u_x_j) ^ 2 +
        (y_j + u_y_i - y_i - u_y_j) ^ 2 +
        (z_j + u_z_i - z_i - u_z_j) ^ 2)

    # Update the element orientation angles at its nodes (i) and (j):
    element.state.ω_i += δu_i[4]
    element.state.ω_j += δu_j[4]

    # Compute the updated element global-to-local transformation matrix:
    Γ′ = compute_Γ(
        x_i + u_x_i, y_i + u_y_i, z_i + u_z_i, element.state.ω_i,
        x_j + u_x_j, y_j + u_y_j, z_j + u_z_j, element.state.ω_j,
        element.state.L)

    # Construct the global element displacement increment vector:
    δu_g = [δu_i; δu_j]

    # Transform the global element displacement increment vector into its previous local coordinate system:
    δu_l = element.state.Γ * δu_g

    # Compute the element internal force increment vector in its previous local coordinate system:
    δq_l  = (element.state.k_e_l + element.state.k_g_l) * δu_l

    # Compute the element internal force increment vector in its global coordinate system:
    q_l   = element.state.q + δq_l
    q_g   = transpose(element.state.Γ) * q_l
    q_l′  = Γ′ * q_g
    element.state.q = q_l′

    element.state.Γ = Γ′

    k_e_l = zeros(T, 12, 12)
    compute_k_e_l!(k_e_l, element, element.state.L)
    condense!(k_e_l, element.releases_i, element.releases_j)
    element.state.k_e_l = k_e_l
    element.state.k_e_g = transform(k_e_l, Γ′)

    k_g_l = zeros(T, 12, 12)
    compute_k_g_l!(k_g_l, element, element.state.L, element.state.q[7])
    condense!(k_g_l, element.releases_i, element.releases_j)
    element.state.k_g_l = k_g_l
    element.state.k_g_g = transform(k_g_l, Γ′)

    element.state.modified = true

    # Return the updated element:
    return element
end

gettype(::Element{NIT, NJT, ST, MT, OAT, ET}) where {NIT, NJT, ST, MT, OAT, ET} = ET

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
