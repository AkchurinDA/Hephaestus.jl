using Test
using Hephaestus
using DifferentiationInterface, ForwardDiff, ReverseDiff

@testset verbose = true "Basics"                     include("Basics.jl")
@testset verbose = true "Automatic differentiation"  include("AutomaticDifferentiation.jl")
@testset verbose = true "Denavit and Hajjar (2013)"  include("Denavit and Hajjar (2013).jl")
@testset verbose = true "Ziemian and Ziemian (2021)" include("Ziemian and Ziemian (2021).jl")
