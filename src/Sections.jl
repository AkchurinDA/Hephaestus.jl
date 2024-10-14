abstract type AbstractSection end

"""
    struct Section

A type representing a section in a finite element model.

$(FIELDS)
"""
struct Section{SPT<:Real} <: AbstractSection
    "Cross-sectional area"
    A   ::SPT
    "Moment of inertia about the local ``z``-axis"
    I_zz::SPT
    "Moment of inertia about the local ``y``-axis"
    I_yy::SPT
    "Polar moment of inertia"
    J   ::SPT
end

Section(A::Real, I_zz::Real, I_yy::Real, J::Real) = Section(promote(A, I_zz, I_yy, J)...)