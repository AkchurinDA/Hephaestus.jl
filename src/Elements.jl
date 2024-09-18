"""
    struct Element

$(FIELDS)
"""
struct Element
    # -------------------------
    # PROVIDED INFORMATION 
    # -------------------------
    "Element tag assigned by the user"
    tag             ::Int
    "Tag of the node (i) of the element"
    node_i_tag      ::Int
    "Tag of the node (j) of the element"
    node_j_tag      ::Int
    "Tag of the material of the element"
    material_tag    ::Int
    "Tag of the section of the element"
    section_tag     ::Int
    "Angle that defines the orientation of the element"
    ω               ::Real

    # ------------------------- 
    # COMPUTED INFORMATION 
    # -------------------------
    "Length of the element"
    L           ::Real
    "Local-to-global transformation matrix of the element"
    T           ::Matrix{<:Real}
    "Element elastic stiffness matrix in the LCS"
    k_e_l       ::Matrix{<:Real}
    "Element elastic stiffness matrix in the GCS"
    k_e_g       ::Matrix{<:Real}
    "Element geometric stiffness matrix in the LCS (without axial force term ``P``)"
    k_g_l       ::Matrix{<:Real}
    "Element geometric stiffness matrix in the GCS (without axial force term ``P``)"
    k_g_g       ::Matrix{<:Real}
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
    ω::EOT,
    L::ELT) where {CTI<:Real, CTJ<:Real, EOT<:Real, ELT<:Real}
    # Compute the rotation angles:
    ρ = -atan(z_j - z_i, x_j - x_i)
    χ = π / 2 - acos((y_j - y_i) / L)

    # Construct the sub-rotation matrix:
    s_ρ, c_ρ = sincos(ρ)
    s_χ, c_χ = sincos(χ)
    s_ω, c_ω = sincos(ω)
    γ = [
        +c_χ * c_ρ                      +s_χ          -c_χ * s_ρ                  ;
        +s_ω * s_ρ - c_ω * s_χ * c_ρ    +c_ω * c_χ    +c_ω * s_χ * s_ρ + s_ω * c_ρ;
        +c_ω * s_ρ + s_ω * s_χ * c_ρ    -s_ω * c_χ    -s_ω * s_χ * s_ρ + c_ω * c_ρ]

    # Remove small values if any:
    map!(x -> abs(x) < eps() ? 0 : x, γ, γ)

    # Preallocate the transformation matrix and fill it:
    T = zeros(eltype(γ), 12, 12)
    @inbounds T[1:3  , 1:3  ] = γ
    @inbounds T[4:6  , 4:6  ] = γ
    @inbounds T[7:9  , 7:9  ] = γ
    @inbounds T[10:12, 10:12] = γ

    # Return the transformation matrix:
    return T
end

#=
using BenchmarkTools
@benchmark _compute_T(
    0.0   , 0.0, 0.0,
    5000.0, 0.0, 0.0,
    π / 4,
    5000.0)

BenchmarkTools.Trial: 10000 samples with 849 evaluations.
    Range (min … max):  130.840 ns …  48.280 μs  ┊ GC (min … max): 0.00% … 99.61%
    Time  (median):     152.925 ns               ┊ GC (median):    0.00%
    Time  (mean ± σ):   165.996 ns ± 482.508 ns  ┊ GC (mean ± σ):  6.10% ±  8.82%
    Memory estimate: 1.48 KiB, allocs estimate: 11.
=#

function _compute_k_e_l(
    E::MPT, ν::MPT,
    A::SPT, I_zz::SPT, I_yy::SPT, J::SPT,
    L::ELT)::Matrix where {MPT<:Real, SPT<:Real, ELT<:Real}
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
    map!(x -> abs(x) < eps() ? 0 : x, k_e_l, k_e_l)
    
    # Return the element stiffness matrix:
    return k_e_l
end

#=
using BenchmarkTools
@benchmark _compute_k_e_l(
    29000.0, 0.3, 
    10.0, 100.0, 100.0, 5.0,
    5000.0)

BenchmarkTools.Trial: 10000 samples with 970 evaluations.
    Range (min … max):  83.977 ns …  42.370 μs  ┊ GC (min … max): 0.00% … 99.59%
    Time  (median):     89.820 ns               ┊ GC (median):    0.00%
    Time  (mean ± σ):   98.868 ns ± 423.162 ns  ┊ GC (mean ± σ):  6.48% ±  7.01%    
    Memory estimate: 1.22 KiB, allocs estimate: 1.
=#

function _compute_k_g_l(
    A::SPT, I_zz::SPT, I_yy::SPT,
    L::ELT) where {SPT<:Real, ELT<:Real}
    # Preallocate:
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
    map!(x -> abs(x) < eps() ? 0 : x, k_g_l, k_g_l)

    # Return the element geometric stiffness matrix:
    return k_g_l
end

#=
using BenchmarkTools
@benchmark _compute_k_g_l(
    10.0, 100.0, 100.0,
    5000.0)

BenchmarkTools.Trial: 10000 samples with 970 evaluations.
    Range (min … max):  83.419 ns …  41.634 μs  ┊ GC (min … max): 0.00% … 99.63%
    Time  (median):     88.360 ns               ┊ GC (median):    0.00%
    Time  (mean ± σ):   97.372 ns ± 415.815 ns  ┊ GC (mean ± σ):  6.49% ±  7.05%
    Memory estimate: 1.22 KiB, allocs estimate: 1.
=#