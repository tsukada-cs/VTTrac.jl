push!(LOAD_PATH,"../src/")

using VTTrac
using Documenter

DocMeta.setdocmeta!(VTTrac, :DocTestSetup, :(using VTTrac); recursive=true)

makedocs(
    sitename = "VTTrac.jl",
    format = Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tsukada-cs.github.io/VTTrac.jl",
    ),
    pages = [
        "Home" => "index.md",
        "How to Use" => "howtouse.md",
        "References" => "references.md",
    ],
)

deploydocs(
    repo = "github.com/tsukada-cs/VTTrac.jl.git",
)
