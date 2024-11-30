struct Section{T<:Real} <: AbstractSection
    A::T
    I_zz::T
    I_yy::T
    J::T

    function Section(A::Real, I_zz::Real, I_yy::Real, J::Real)
        T = float(promote_type(
            typeof(A), 
            typeof(I_zz), 
            typeof(I_yy), 
            typeof(J)))

        return new{T}(A, I_zz, I_yy, J)
    end
end