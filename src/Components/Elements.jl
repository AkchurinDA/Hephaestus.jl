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
    ID        ::Int
    node_i    ::Node{NIT}
    node_j    ::Node{NJT}
    section   ::Section{ST}
    material  ::Material{MT}
    ω         ::OAT
    releases_i::Vector{Bool}
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

mutable struct ElementState
    ID   ::Int
    x_i  ::Real
    y_i  ::Real
    z_i  ::Real
    x_j  ::Real
    y_j  ::Real
    z_j  ::Real
    L    ::Real
    ω_i  ::Real
    ω_j  ::Real
    Γ    ::AbstractMatrix{<:Real}
    q    ::AbstractVector{<:Real}
    k_e_l::AbstractMatrix{<:Real}
    k_e_g::AbstractMatrix{<:Real}
    k_g_l::AbstractMatrix{<:Real}
    k_g_g::AbstractMatrix{<:Real}
end

function compute_γ(β::Real, χ::Real, ω::Real)
    s_β, c_β = sincos(β)
    s_χ, c_χ = sincos(χ)
    s_ω, c_ω = sincos(ω)

    γ = [
        +c_χ * c_β                      +s_χ          -c_χ * s_β                  ;
        -c_ω * s_χ * c_β + s_ω * s_β    +c_ω * c_χ    +c_ω * s_χ * s_β + s_ω * c_β;
        +s_ω * s_χ * c_β + c_ω * s_β    -s_ω * c_χ    -s_ω * s_χ * s_β + c_ω * c_β]

    return γ
end

@memoize function compute_Γ(
    x_i::Real, y_i::Real, z_i::Real,
    x_j::Real, y_j::Real, z_j::Real,
    ω_i::Real,
    ω_j::Real)::AbstractMatrix{<:Real}
    # Compute the element length projections:
    Δx = x_j - x_i
    Δy = y_j - y_i
    Δz = z_j - z_i

    # Compute the element length:
    L = sqrt(Δx ^ 2 + Δy ^ 2 + Δz ^ 2)

    # Compute the element orientation angles:
    β = -atan(Δz, Δx) # Conventionally, this angle is called "ρ", but had to rename it to avoid conflicts with the "ρ" symbol used for the material's density
    χ = π / 2 - acos(Δy / L)

    # Compute the element subtransformation matrix:
    γ_i = compute_γ(β, χ, ω_i)
    γ_j = compute_γ(β, χ, ω_j)

    # Compute the element transformation matrix:
    T = promote_type(eltype(γ_i), eltype(γ_j))
    Γ = zeros(T, 12, 12)
    Γ[1:3  , 1:3  ] .= γ_i
    Γ[4:6  , 4:6  ] .= γ_i
    Γ[7:9  , 7:9  ] .= γ_j
    Γ[10:12, 10:12] .= γ_j

    map!(x -> abs(x) ≤ eps() ? zero(T) : x, Γ, Γ)

    # Return the element transformation matrix:
    return Γ
end

function compute_k_e_l!(k_e_l::AbstractMatrix{<:Real}, element::Element, L::Real)::AbstractMatrix{<:Real}
    # Extract the section properties:
    E, ν = element.material.E, element.material.ν

    # Extract the section properties:
    A, I_zz, I_yy, J = element.section.A, element.section.I_zz, element.section.I_yy, element.section.J

    # Compute the element elastic stiffness matrix:
    compute_k_e_l!(k_e_l, E, ν, A, I_zz, I_yy, J, L)

    # Return the element elastic stiffness matrix:
    return k_e_l
end

@memoize function compute_k_e_l!(k_e_l::AbstractMatrix{<:Real},
    E::Real, ν::Real, 
    A::Real, I_zz::Real, I_yy::Real, J::Real, 
    L::Real)::AbstractMatrix{<:Real}
    # Compute the shear modulus:
    G = E / (2 * (1 + ν))

    # Compute the element elastic stiffness matrix:
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

    for i in 1:12, j in (i + 1):12
        @inbounds k_e_l[j, i] = k_e_l[i, j]
    end

    map!(x -> abs(x) ≤ eps() ? zero(eltype(k_e_l)) : x, k_e_l, k_e_l)

    # Return the element elastic stiffness matrix:
    return k_e_l
end

function compute_k_g_l!(k_g_l::AbstractMatrix{<:Real}, element::Element, L::Real, N::Real)::AbstractMatrix{<:Real}
    # Extract the section properties:
    A, J = element.section.A, element.section.J

    # Compute the element elastic stiffness matrix:
    compute_k_g_l!(k_g_l, A, J, L, N)

    # Return the element elastic stiffness matrix:
    return k_g_l
end

@memoize function compute_k_g_l!(k_g_l::AbstractMatrix{<:Real},
    A::Real, J::Real,
    L::Real,
    N::Real)::AbstractMatrix{<:Real}
    # Compute the element geometric stiffness matrix:
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

    for i in 1:12, j in (i + 1):12
        @inbounds k_g_l[j, i] = k_g_l[i, j]
    end

    k_g_l .*= N

    map!(x -> abs(x) ≤ eps() ? zero(eltype(k_g_l)) : x, k_g_l, k_g_l)

    # Return the element geometric stiffness matrix:
    return k_g_l
end

@memoize function compute_m_l!(m_l::AbstractMatrix{<:Real},
    ρ::Real, 
    A::Real, J::Real, 
    L::Real)::AbstractMatrix{<:Real}
    # Compute the element mass matrix:
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
    
    for i in 1:12, j in (i + 1):12
        @inbounds m_l[j, i] = m_l[i, j]
    end

    m_l .*= (ρ * A * L) / 420

    # Return the element mass matrix:
    return m_l
end

@memoize function condense!(m::AbstractMatrix{<:Real}, releases_i::Vector{Bool}, releases_j::Vector{Bool})::AbstractMatrix{<:Real}
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

@memoize function transform(m::AbstractMatrix{<:Real}, Γ::AbstractMatrix{<:Real})::AbstractMatrix{<:Real}
    # Transform the matrix to the global coordinate system:
    M = transpose(Γ) * m * Γ

    # Return the transformed matrix:
    return M
end

gettype(::Element{NIT, NJT, ST, MT, OAT, ET}) where {NIT, NJT, ST, MT, OAT, ET} = ET