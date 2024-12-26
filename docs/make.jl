push!(LOAD_PATH, "../src/")

using Hephaestus
using Documenter

makedocs(
    sitename = "Hephaestus.jl",
    authors  = "Damir Akchurin",
    pages    = [
        "Home" => "index.md",
        "API"  => "API.md"])

deploydocs(
    repo = "github.com/AkchurinDA/Hephaestus.jl")