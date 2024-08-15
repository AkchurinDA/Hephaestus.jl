struct Section{T<:Real}
    "Section ID."
    ID::Integer
    "Gross cross-sectional area."
    A::T
    "Moment of inertia about the local z-axis."
    I::T

    function Section(ID::Integer, A::Real, I::Real)
        # Promote the section properties to the same common type:
        section_properties = promote(A, I)

        # Construct the section:
        return new{eltype(section_properties)}(ID, section_properties...)
    end
end