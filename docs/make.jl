push!(LOAD_PATH, "../src/")

using Hephaestus
using Documenter

makedocs(
    sitename = "Hephaestus.jl",
    authors  = "Damir Akchurin",
    pages    = [
        "Home"        => "index.md",
        "Quick Start" => "QuickStart.md",
        "Documentation" => [
            "Section library"  => "SectionLibrary.md",
            "Material library" => "MaterialLibrary.md",
            "Element library"  => "ElementLibrary.md"],
        "API"         => "API.md"])

deploydocs(
    repo = "github.com/AkchurinDA/Hephaestus.jl")