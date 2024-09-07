using Hephaestus
using Documenter
using DocumenterCitations

makedocs(
    sitename = "Hephaestus.jl",
    authors = "Damir Akchurin, AkchurinDA@gmail.com",
    pages = [
        "Home" => "index.md"
        "Quick Start" => "QuickStart.md"],
    plugins = [
        CitationBibliography(
            joinpath(@__DIR__, "src", "References.bib"),
            style=:numeric)])

deploydocs(
    repo = "github.com/AkchurinDA/Hephaestus.jl")