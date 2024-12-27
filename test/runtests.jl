using Test
using Hephaestus
using DifferentiationInterface, ForwardDiff, ReverseDiff

@testset "Denavit and Hajjar (2013).jl" begin
    include("Denavit and Hajjar (2013).jl")
end

@testset "Ziemian and Ziemian (2021).jl" begin
    include("Ziemian and Ziemian (2021).jl")
end