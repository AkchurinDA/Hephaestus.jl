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
    L         ::Real
    Γ         ::AbstractMatrix{ET}
    k_e_l     ::AbstractMatrix{ET}
    k_e_g     ::AbstractMatrix{ET}
    k_g_l     ::AbstractMatrix{ET}
    k_g_g     ::AbstractMatrix{ET}
    m_l       ::AbstractMatrix{ET}
    m_g       ::AbstractMatrix{ET}

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

        # Extract the coordinates of the nodes:
        x_i, y_i, z_i = node_i.x, node_i.y, node_i.z
        x_j, y_j, z_j = node_j.x, node_j.y, node_j.z

        # Compute the element length projections:
        Δx = x_j - x_i
        Δy = y_j - y_i
        Δz = z_j - z_i

        # Compute the element length:
        L = sqrt(Δx ^ 2 + Δy ^ 2 + Δz ^ 2)

        # Compute the transformation matrix:
        Γ = compute_Γ(x_i, y_i, z_i, x_j, y_j, z_j, ω)

        # Extract the section properties:
        A, I_zz, I_yy, J = section.A, section.I_zz, section.I_yy, section.J

        # Extract the material properties:
        E, ν, ρ = material.E, material.ν, material.ρ

        # Compute the element elastic stiffness matrix in the local coordinate system:
        k_e_l = zeros(ET, 12, 12)
        compute_k_e_l!(k_e_l, E, ν, A, I_zz, I_yy, J, L)
        condense!(k_e_l, releases_i, releases_j)
        k_e_g = transform(k_e_l, Γ)

        # Compute the element geometric stiffness matrix in the local coordinate system:
        k_g_l = zeros(ET, 12, 12)
        compute_k_g_l!(k_g_l, A, J, L)
        condense!(k_g_l, releases_i, releases_j)
        k_g_g = transform(k_g_l, Γ)

        # Compute the element mass matrix in the local coordinate system:
        m_l = zeros(ET, 12, 12)
        compute_m_l!(m_l, ρ, A, J, L)
        condense!(m_l, releases_i, releases_j)
        m_g = transform(m_l, Γ)

        # Return the element:
        return new{NIT, NJT, ST, MT, OAT, ET}(ID, node_i, node_j, section, material, ω, releases_i, releases_j, L, Γ, k_e_l, k_e_g, k_g_l, k_g_g, m_l, m_g)
    end
end

mutable struct ElementState
    ID   ::Int
    q    ::AbstractVector{<:Real}
    x_i  ::Real
    y_i  ::Real
    z_i  ::Real
    x_j  ::Real
    y_j  ::Real
    z_j  ::Real
    ω_i  ::Real
    ω_j  ::Real
    Γ    ::AbstractMatrix{<:Real}
    k_e_l::AbstractMatrix{<:Real}
    k_e_g::AbstractMatrix{<:Real}
    k_g_l::AbstractMatrix{<:Real}
    k_g_g::AbstractMatrix{<:Real}
end

# TODO: Add a function to update the element state based of the current global displacement vector increment.

@memoize function compute_Γ(
    x_i::Real, y_i::Real, z_i::Real,
    x_j::Real, y_j::Real, z_j::Real,
    ω::Real)::AbstractMatrix{<:Real}
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
    s_β, c_β = sincos(β)
    s_χ, c_χ = sincos(χ)
    s_ω, c_ω = sincos(ω)
    γ = [
        +c_χ * c_β                      +s_χ          -c_χ * s_β                  ;
        -c_ω * s_χ * c_β + s_ω * s_β    +c_ω * c_χ    +c_ω * s_χ * s_β + s_ω * c_β;
        +s_ω * s_χ * c_β + c_ω * s_β    -s_ω * c_χ    -s_ω * s_χ * s_β + c_ω * c_β]

    # Compute the element transformation matrix:
    Γ = zeros(eltype(γ), 12, 12)
    Γ[1:3  , 1:3  ] .= γ
    Γ[4:6  , 4:6  ] .= γ
    Γ[7:9  , 7:9  ] .= γ
    Γ[10:12, 10:12] .= γ

    # Return the element transformation matrix:
    return Γ
end

function compute_Γ(
    x_i::Real, y_i::Real, z_i::Real,
    x_j::Real, y_j::Real, z_j::Real,
    u_x_i::Real, u_y_i::Real, u_z_i::Real, θ_x_i::Real, θ_y_i::Real, θ_z_i::Real,
    u_x_j::Real, u_y_j::Real, u_z_j::Real, θ_x_j::Real, θ_y_j::Real, θ_z_j::Real,
    ω_i::Real,
    ω_j::Real)::AbstractMatrix{<:Real}
    # Compute the element length projections:
    Δx = (x_j + u_x_j) - (x_i + u_x_i)
    Δy = (y_j + u_y_j) - (y_i + u_y_i)
    Δz = (z_j + u_z_j) - (z_i + u_z_i)

    # Compute the element length:
    L = sqrt(Δx ^ 2 + Δy ^ 2 + Δz ^ 2)

    # Compute the element orientation angles:
    β = -atan(Δz, Δx) # Conventionally, this angle is called "ρ", but had to rename it to avoid conflicts with the "ρ" symbol used for the material's density
    χ = π / 2 - acos(Δy / L)

    # Compute the element subtransformation matrix for node (i):
    β_i = β + θ_y_i
    χ_i = χ + θ_z_i
    ω_i = ω_i + θ_x_i
    s_β_i, c_β_i = sincos(β_i)
    s_χ_i, c_χ_i = sincos(χ_i)
    s_ω_i, c_ω_i = sincos(ω_i)
    γ_i = [
        +c_χ_i * c_β_i                            +s_χ_i            -c_χ_i * s_β_i                        ;
        -c_ω_i * s_χ_i * c_β_i + s_ω_i * s_β_i    +c_ω_i * c_χ_i    +c_ω_i * s_χ_i * s_β_i + s_ω_i * c_β_i;
        +s_ω_i * s_χ_i * c_β_i + c_ω_i * s_β_i    -s_ω_i * c_χ_i    -s_ω_i * s_χ_i * s_β_i + c_ω_i * c_β_i]

    # Compute the element subtransformation matrix for node (j):
    β_j = β + θ_y_j
    χ_j = χ + θ_z_j
    ω_j = ω_j + θ_x_j
    s_β_j, c_β_j = sincos(β_j)
    s_χ_j, c_χ_j = sincos(χ_j)
    s_ω_j, c_ω_j = sincos(ω_j)
    γ_j = [
        +c_χ_j * c_β_j                            +s_χ_j            -c_χ_j * s_β_j                        ;
        -c_ω_j * s_χ_j * c_β_j + s_ω_j * s_β_j    +c_ω_j * c_χ_j    +c_ω_j * s_χ_j * s_β_j + s_ω_j * c_β_j;
        +s_ω_j * s_χ_j * c_β_j + c_ω_j * s_β_j    -s_ω_j * c_χ_j    -s_ω_j * s_χ_j * s_β_j + c_ω_j * c_β_j]

    # Compute the element transformation matrix:
    T = promote_type(eltype(γ_i), eltype(γ_j))
    Γ = zeros(T, 12, 12)
    Γ[1:3  , 1:3  ] .= γ_i
    Γ[4:6  , 4:6  ] .= γ_i
    Γ[7:9  , 7:9  ] .= γ_j
    Γ[10:12, 10:12] .= γ_j

    # Return the element transformation matrix:
    return Γ
end

@memoize function compute_k_e_l!(k_e_l::AbstractMatrix{T}, # element::Element)::AbstractMatrix{T} where {T <: Real}
    E::MT, ν::MT, 
    A::ST, I_zz::ST, I_yy::ST, J::ST, 
    L::LT)::AbstractMatrix{T} where {
        T  <: Real,
        MT <: Real, 
        ST <: Real, 
        LT <: Real}
    # # Extract the material properties:
    # E, ν = element.material.E, element.material.ν

    # # Extract the section properties:
    # A, I_zz, I_yy, J = element.section.A, element.section.I_zz, element.section.I_yy, element.section.J

    # # Extract the element length:
    # L = element.L

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

    # Return the element elastic stiffness matrix:
    return k_e_l
end

@memoize function compute_k_g_l!(k_g_l::AbstractMatrix{T}, # element::Element)::AbstractMatrix{T} where {T <: Real}
    A::ST, J::ST,
    L::LT)::AbstractMatrix{T} where {
        T  <: Real,
        ST <: Real, 
        LT <: Real}
    # # Extract the section properties:
    # A, J = element.section.A, element.section.J

    # # Extract the element length:
    # L = element.L

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

    # Return the element geometric stiffness matrix:
    return k_g_l
end

@memoize function compute_m_l!(m_l::AbstractMatrix{T},
    ρ::MT, 
    A::ST, J::ST, 
    L::LT)::AbstractMatrix{T} where {
        T  <: Real,
        MT <: Real, 
        ST <: Real, 
        LT <: Real}
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

@memoize function transform(m::AbstractMatrix{<:Real}, Γ::AbstractMatrix{<:Real})
    # Transform the matrix to the global coordinate system:
    M = Γ' * m * Γ

    # Return the transformed matrix:
    return M
end

gettype(::Element{NIT, NJT, ST, MT, OAT, ET}) where {NIT, NJT, ST, MT, OAT, ET} = ET