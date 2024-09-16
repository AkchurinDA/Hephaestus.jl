"""
    struct Section

A type that represents a section in the FE model of a structure.

This type should not be used directly. 
Instead, use the [`add_section!()`](@ref) function to add a section to the model.
To remove a section from the model, use the [`del_section!()`](@ref) function.
"""
struct Section{SPT <: Real}
    ID_u    ::Int
    A       ::SPT
    I_zz    ::SPT
    I_yy    ::SPT
    J       ::SPT
end

Section(ID::Int, A::Real, I_zz::Real, I_yy::Real, J::Real) = Section(ID, promote(A, I_zz, I_yy, J)...)