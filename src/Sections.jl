"""
    Section

A type that represents a section in the finite element model of a structure.

# Fields
$(FIELDS)
"""
struct Section{SPT<:Real}
    "Unique identifier of the section provided by the user"
    ID      ::Int
    "Cross-sectional area of the section, ``A``"
    A       ::SPT
    "Moment of inertia about the local ``z``-axis of the section, ``I_{zz}``"
    I_zz    ::SPT
    "Moment of inertia about the local ``y``-axis of the section, ``I_{yy}``"
    I_yy    ::SPT
    "Polar moment of inertia of the section, ``J``"
    J       ::SPT
end
Section(ID::Int, A::Real, I_zz::Real, I_yy::Real, J::Real) = Section(ID, promote(A, I_zz, I_yy, J)...)