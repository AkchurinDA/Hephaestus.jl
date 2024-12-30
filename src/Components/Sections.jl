"""
    struct Section

A type representing a section in a finite element model.

$(FIELDS)
"""
struct Section{T<:Real}
    "Identification tag"
    ID::Int
    "Cross-sectional area, ``A``"
    A::T
    "Moment of inertia about the local ``z``-axis, ``I_{zz}``"
    I_zz::T
    "Moment of inertia about the local ``y``-axis, ``I_{yy}``"
    I_yy::T
    "Torsional constant, ``J``"
    J::T

    function Section(ID::Int, 
        A::T1, I_zz::T2, I_yy::T3, J::T4) where {
        T1<:Real,
        T2<:Real,
        T3<:Real,
        T4<:Real}
        # Promote the types:
        T = float(promote_type(T1, T2, T3, T4))

        # Construct a new instance:
        return new{T}(ID, A, I_zz, I_yy, J)
    end
end

gettype(::Section{T}) where {T} = T