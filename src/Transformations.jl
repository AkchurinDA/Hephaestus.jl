function compute_T(
    x_i::NICT, y_i::NICT, z_i::NICT,
    x_j::NJCT, y_j::NJCT, z_j::NJCT,
    ω::EOT,
    L::EPT) where {NICT,NJCT,EOT,EPT}
    # Preallocate the transformation matrix:
    T = SparseArrays.spzeros(float(promote_type(NICT, NJCT, EOT, EPT)), 12, 12)

    # Construct the sub-transformation matrix:
    ρ = -atan(z_j - z_i, x_j - x_i)
    χ = π / 2 - acos((y_j - y_i) / L)
    γ_ρ = [+cos(ρ) 0 -sin(ρ); 0 +1 0; +sin(ρ) 0 +cos(ρ)]
    γ_χ = [+cos(χ) +sin(χ) 0; -sin(χ) +cos(χ) 0; 0 0 +1]
    γ_ω = [+1 0 0 ;0 +cos(ω) +sin(ω); 0 -sin(ω) +cos(ω)]
    γ = γ_ω * γ_χ * γ_ρ

    # Define the elements of the transformation matrix:
    @inline T[  1:3,   1:3] .= γ
    @inline T[  4:6,   4:6] .= γ
    @inline T[  7:9,   7:9] .= γ
    @inline T[10:12, 10:12] .= γ

    # Return the result:
    return T
end