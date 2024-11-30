struct CLoad{T<:Real}
    F_x::T
    F_y::T
    F_z::T
    M_x::T
    M_y::T
    M_z::T

    function CLoad(F_x::Real, F_y::Real, F_z::Real, M_x::Real, M_y::Real, M_z::Real)
        T = float(promote_type(
            typeof(F_x), 
            typeof(F_y), 
            typeof(F_z), 
            typeof(M_x), 
            typeof(M_y), 
            typeof(M_z)))

        return new{T}(F_x, F_y, F_z, M_x, M_y, M_z)
    end
end

struct DLoad{T<:Real}
    q_x::T
    q_y::T
    q_z::T

    function DLoad(q_x::Real, q_y::Real, q_z::Real)
        T = float(promote_type(
            typeof(q_x), 
            typeof(q_y), 
            typeof(q_z)))

        return new{T}(q_x, q_y, q_z)
    end
end