"""
    struct Section

A type representing a section in the model of a structure of interest.

$(FIELDS)
"""
struct Section{T<:Real}
    "Unique identifier"
    ID      ::Int
    "Cross-sectional area, ``A``"
    A       ::T
    "Moment of inertia about the local ``z``-axis, ``I_{zz}``"
    I_zz    ::T
    "Moment of inertia about the local ``y``-axis, ``I_{yy}``"
    I_yy    ::T
    "Polar moment of inertia about the local ``x``-axis, ``J``"
    J       ::T
end

Section(ID, A, I_zz, I_yy, J) = Section(ID, promote(A, I_zz, I_yy, J)...)