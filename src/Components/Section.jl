struct Section{T<:Real}
    "Section ID."
    ID::Integer
    "Gross cross-sectional area."
    A::T
    "Moment of inertia about the z-axis."
    I_zz::T
    "Moment of inertia about the y-axis."
    I_yy::T
    "Polar moment of inertia."
    J::T

    function Section(ID::Integer, A::Real, I_zz::Real, I_yy::Real, J::Real)
        # Promote the section properties to the same common type:
        SectionProperties = promote(A, I_zz, I_yy, J)

        # Construct the section:
        return new{eltype(SectionProperties)}(ID, SectionProperties...)
    end
end