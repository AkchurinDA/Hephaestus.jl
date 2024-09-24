"""
    struct Section

A type representing a section in the finite element model of a structure of interest.
"""
struct Section{T<:Real}
    ID      ::Int
    A       ::T
    I_zz    ::T
    I_yy    ::T
    J       ::T
end

Section(ID::Int, A::Real, I_zz::Real, I_yy::Real, J::Real) = Section(ID, promote(A, I_zz, I_yy, J)...)