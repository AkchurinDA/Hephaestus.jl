"""
    struct Model

A type representing the finite element model of a structure of interest.
"""
@kwdef struct Model
    nodes       ::Dict{Int, Node    } = Dict{Int, Node    }()
    materials   ::Dict{Int, Material} = Dict{Int, Material}()
    sections    ::Dict{Int, Section } = Dict{Int, Section }()
    elements    ::Dict{Int, Element } = Dict{Int, Element }()

    supports    ::Dict{Int, Vector{Bool}} = Dict{Int, Vector{Bool}}()

    conc_loads  ::Dict{Int, Vector{<:Real}} = Dict{Int, Vector{<:Real}}()
    dist_loads  ::Dict{Int, Vector{<:Real}} = Dict{Int, Vector{<:Real}}()

    p_l         ::Dict{Int, Vector{<:Real}} = Dict{Int, Vector{<:Real}}()
    p_g         ::Dict{Int, Vector{<:Real}} = Dict{Int, Vector{<:Real}}()
end