struct Section{SPT<:Real}
    tag     ::Int
    A       ::SPT
    I_zz    ::SPT
    I_yy    ::SPT
    J       ::SPT
end

Section(tag::Int, A::Real, I_zz::Real, I_yy::Real, J::Real) = Section(tag, promote(A, I_zz, I_yy, J)...)