using VTTrac
using Documenter

DocMeta.setdocmeta!(VTTrac, :DocTestSetup, :(using VTTrac); recursive=true)

makedocs(;
    modules=[VTTrac],
    authors="tsukada-cs <tsukada.cs@gmail.com> and contributors",
    repo="https://github.com/tsukada-cs/VTTrac.jl/blob/{commit}{path}#{line}",
    sitename="VTTrac.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tsukada-cs.github.io/VTTrac.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tsukada-cs/VTTrac.jl",
    devbranch="main",
)
