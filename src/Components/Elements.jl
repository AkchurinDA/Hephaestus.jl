struct Element{NIT<:Real, NJT<:Real, ST<:Real, MT<:Real}
    "Identification tag"
    ID::Int
    "Node (``i``)"
    node_i::Node{NIT}
    "Node (``j``)"
    node_j::Node{NJT}
    "Section"
    section::Section{ST}
    "Material"
    material::Material{MT}
    "Orientation angle, ``\\omega``"
    ω::Real
    releases_i::Vector{Bool}
    releases_j::Vector{Bool}
    L::Real
    Γ::Matrix{<:Real}
    k_e_l::Matrix{<:Real}
    k_g_l::Matrix{<:Real}
    m_l::Matrix{<:Real}

    function Element(ID, 
        node_i::Node{NIT}, node_j::Node{NJT},
        section::Section{ST}, 
        material::Material{MT},
        ω::Real,
        releases_i::Vector{<:Bool},
        releases_j::Vector{<:Bool}) where {
        NIT<:Real,
        NJT<:Real,
        ST<:Real,
        MT<:Real}
        # Extract the coordinates of the nodes:
        x_i, y_i, z_i = node_i.x, node_i.y, node_i.z
        x_j, y_j, z_j = node_j.x, node_j.y, node_j.z

        # Extract the section properties:
        A, I_zz, I_yy, J = section.A, section.I_zz, section.I_yy, section.J

        # Extract the material properties:
        E, ν, ρ = material.E, material.ν, material.ρ

        # Compute the element length projections:
        Δx = x_j - x_i
        Δy = y_j - y_i
        Δz = z_j - z_i

        # Compute the element length:
        L = sqrt(Δx ^ 2 + Δy ^ 2 + Δz ^ 2)

        # Compute the element orientation angles:
        ρ = -atan(Δz / Δx)
        χ = π / 2 - acos(Δy / L)

        # Compute the element subtransformation matrix:
        s_ρ, c_ρ = sincos(ρ)
        s_χ, c_χ = sincos(χ)
        s_ω, c_ω = sincos(ω)
        γ = [
            +c_χ * c_ρ                      +s_χ          -c_χ * s_ρ                  ;
            -c_ω * s_χ * c_ρ + s_ω * s_ρ    +c_ω * c_χ    +c_ω * s_χ * s_ρ + s_ω * c_ρ;
            +s_ω * s_χ * c_ρ + c_ω * s_ρ    -s_ω * c_χ    -s_ω * s_χ * s_ρ + c_ω * c_ρ]

        # Compute the element transformation matrix:
        Γ = zeros(eltype(γ), 12, 12)
        Γ[1:3  , 1:3  ] = γ
        Γ[4:6  , 4:6  ] = γ
        Γ[7:9  , 7:9  ] = γ
        Γ[10:12, 10:12] = γ

        # Compute the element elastic stiffness matrix in the local coordinate system:
        k_e_l = compute_k_e_l(E, ν, A, I_zz, I_yy, J, L)
        condense!(k_e_l, releases_i, releases_j)

        # Compute the element geometric stiffness matrix in the local coordinate system:
        k_g_l = compute_k_g_l(A, J, L)
        condense!(k_g_l, releases_i, releases_j)

        # Compute the element mass matrix in the local coordinate system:
        m_l = compute_m_l(ρ, A, J, L)
        condense!(m_l, releases_i, releases_j)

        # Return the element:
        return new{NIT, NJT, ST, MT}(ID, node_i, node_j, section, material, ω, releases_i, releases_j, L, Γ, k_e_l, k_g_l, m_l)
    end
end

@memoize function compute_k_e_l(
    E::Real, ν::Real, 
    A::Real, I_zz::Real, I_yy::Real, J::Real, 
    L::Real)
    # Initialize the element elastic stiffness matrix:
    k_e_l = zeros(Real, 12, 12)

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

@memoize function compute_k_g_l(
    A::Real, J::Real,
    L::Real)
    # Initialize the element geometric stiffness matrix:
    k_g_l = zeros(Real, 12, 12)

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

@memoize function compute_m_l(
    ρ::Real, 
    A::Real, J::Real, 
    L::Real)
    # Initialize the element mass matrix:
    m_l = zeros(Real, 12, 12)

    # Compute the element mass matrix:
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

    # Return the element mass matrix:
    return m_l
end

@memoize function condense!(m::Matrix{<:Real}, releases_i::Vector{Bool}, releases_j::Vector{Bool})
    # Condense the matrix if end releases are present:
    for i in 1:6
        # Node (i):
        if releases_i[i]
            m[i, :] .= 0
            m[:, i] .= 0
        end

        # Node (j):
        if releases_j[i]
            m[6 + i, :] .= 0
            m[:, 6 + i] .= 0
        end
    end

    # Return the condensed matrix:
    return m
end