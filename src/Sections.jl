abstract type AbstractSection end

struct Section{SPT} <: AbstractSection
    ID
    A::SPT
    I_zz::SPT
    I_yy::SPT
    J::SPT

    function Section(ID, A::Real, I_zz::Real, I_yy::Real, J::Real)
        section_properties = promote(A, I_zz, I_yy, J)

        return new{eltype(section_properties)}(ID, section_properties...)
    end
end