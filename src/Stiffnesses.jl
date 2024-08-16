function compute_k_e_l(
    # Material properties:
    E::MPT, ν::MPT,
    # Section properties:
    A::SPT, I_zz::SPT, I_yy::SPT, J::SPT,
    # Element properties:
    L::EPT) where {MPT,SPT,EPT}
    # Preallocate the elastic stiffness matrix:
    k_e_l = SparseArrays.spzeros(float(promote_type(MPT, SPT, EPT)), 12, 12)

    # Define the upper triangular part of the elastic stiffness matrix:
    @inline k_e_l[ 1,  1] = +E * A / L
    @inline k_e_l[ 1,  7] = -E * A / L
    @inline k_e_l[ 2,  2] = +12 * E * I_zz / L^3
    @inline k_e_l[ 2,  6] = +6 * E * I_zz / L^2
    @inline k_e_l[ 2,  8] = -12 * E * I_zz / L^3
    @inline k_e_l[ 2, 12] = +6 * E * I_zz / L^2
    @inline k_e_l[ 3,  3] = +12 * E * I_yy / L^3
    @inline k_e_l[ 3,  5] = -6 * E * I_yy / L^2
    @inline k_e_l[ 3,  9] = -12 * E * I_yy / L^3
    @inline k_e_l[ 3, 11] = -6 * E * I_yy / L^2
    @inline k_e_l[ 4,  4] = +E * J / (2 * (1 + ν) * L)
    @inline k_e_l[ 4, 10] = -E * J / (2 * (1 + ν) * L)
    @inline k_e_l[ 5,  5] = +4 * E * I_yy / L
    @inline k_e_l[ 5,  9] = +6 * E * I_yy / L^2
    @inline k_e_l[ 5, 11] = +2 * E * I_yy / L
    @inline k_e_l[ 6,  6] = +4 * E * I_zz / L
    @inline k_e_l[ 6,  8] = -6 * E * I_zz / L^2
    @inline k_e_l[ 6, 12] = +2 * E * I_zz / L
    @inline k_e_l[ 7,  7] = +E * A / L
    @inline k_e_l[ 8,  8] = -12 * E * I_zz / L^3
    @inline k_e_l[ 8, 12] = -6 * E * I_zz / L^2
    @inline k_e_l[ 9,  9] = +12 * E * I_yy / L^3
    @inline k_e_l[ 9, 11] = +6 * E * I_yy / L^2
    @inline k_e_l[10, 10] = +E * J / (2 * (1 + ν) * L)
    @inline k_e_l[11, 11] = +12 * E * I_yy / L^3
    @inline k_e_l[12, 12] = +4 * E * I_zz / L

    # Define the lower triangular part of the elastic stiffness matrix:
    @inline k_e_l[ 5, 3] = k_e_l[ 3, 5]
    @inline k_e_l[ 6, 2] = k_e_l[ 2, 6]
    @inline k_e_l[ 7, 1] = k_e_l[ 1, 7]
    @inline k_e_l[ 8, 2] = k_e_l[ 2, 8]
    @inline k_e_l[ 8, 6] = k_e_l[ 6, 8]
    @inline k_e_l[ 9, 3] = k_e_l[ 3, 9]
    @inline k_e_l[ 9, 5] = k_e_l[ 5, 9]
    @inline k_e_l[10, 4] = k_e_l[4, 10]
    @inline k_e_l[11, 3] = k_e_l[3, 11]
    @inline k_e_l[11, 5] = k_e_l[5, 11]
    @inline k_e_l[11, 9] = k_e_l[9, 11]
    @inline k_e_l[12, 2] = k_e_l[2, 12]
    @inline k_e_l[12, 6] = k_e_l[6, 12]
    @inline k_e_l[12, 8] = k_e_l[8, 12]

    # Return the result:
    return k_e_l
end

function compute_k_g_l(
    # Loads:
    P::LT,
    # Section properties:
    A::SPT, I_zz::SPT, I_yy::SPT,
    # Element properties:
    L::EPT) where {LT,SPT,EPT}
    # Preallocate the geometric stiffness matrix:
    k_g_l = SparseArrays.spzeros(float(promote_type(LT, SPT, EPT)), 12, 12)

    # Define the upper triangular part of the geometric stiffness matrix:
    @inline k_g_l[ 1,  1] = +P / L
    @inline k_g_l[ 1,  7] = -P / L
    @inline k_g_l[ 2,  2] = +6 * P / (5 * L)
    @inline k_g_l[ 2,  6] = +P / 10
    @inline k_g_l[ 2,  8] = -6 * P / (5 * L)
    @inline k_g_l[ 2, 12] = +P / 10
    @inline k_g_l[ 3,  3] = +6 * P / (5 * L)
    @inline k_g_l[ 3,  5] = -P / 10
    @inline k_g_l[ 3,  9] = -6 * P / (5 * L)
    @inline k_g_l[ 3, 11] = -P / 10
    @inline k_g_l[ 4,  4] = +P * (I_zz + I_yy) / (A * L)
    @inline k_g_l[ 4, 10] = -P * (I_zz + I_yy) / (A * L)
    @inline k_g_l[ 5,  5] = +2 * P * L / 15
    @inline k_g_l[ 5,  9] = +P / 10
    @inline k_g_l[ 5, 11] = -P * L / 30
    @inline k_g_l[ 6,  6] = +2 * P * L / 15
    @inline k_g_l[ 6,  8] = -P / 10
    @inline k_g_l[ 6, 12] = -P * L / 30
    @inline k_g_l[ 7,  7] = +P / L
    @inline k_g_l[ 8,  8] = +6 * P / (5 * L)
    @inline k_g_l[ 8, 12] = -P / 10
    @inline k_g_l[ 9,  9] = +6 * P / (5 * L)
    @inline k_g_l[ 9, 11] = +P / 10
    @inline k_g_l[10, 10] = +P * (I_zz + I_yy) / (A * L)
    @inline k_g_l[11, 11] = +2 * P * L / 15
    @inline k_g_l[12, 12] = +2 * P * L / 15

    # Define the lower triangular part of the geometric stiffness matrix:
    @inline k_g_l[ 5, 3] = k_g_l[3,  5]
    @inline k_g_l[ 6, 2] = k_g_l[2,  6]
    @inline k_g_l[ 7, 1] = k_g_l[1,  7]
    @inline k_g_l[ 8, 2] = k_g_l[2,  8]
    @inline k_g_l[ 8, 6] = k_g_l[6,  8]
    @inline k_g_l[ 9, 3] = k_g_l[3,  9]
    @inline k_g_l[ 9, 5] = k_g_l[5,  9]
    @inline k_g_l[10, 4] = k_g_l[4, 10]
    @inline k_g_l[11, 3] = k_g_l[3, 11]
    @inline k_g_l[11, 5] = k_g_l[5, 11]
    @inline k_g_l[11, 9] = k_g_l[9, 11]
    @inline k_g_l[12, 2] = k_g_l[2, 12]
    @inline k_g_l[12, 6] = k_g_l[6, 12]
    @inline k_g_l[12, 8] = k_g_l[8, 12]

    # Return the result:
    return k_g_l
end