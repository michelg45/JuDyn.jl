using JuDyn
using Documenter

DocMeta.setdocmeta!(JuDyn, :DocTestSetup, :(using JuDyn); recursive=true)

makedocs(;
    modules=[JuDyn],
    authors="Michel Geradin <mgeradin@gmail.com> and contributors",
    repo="https://github.com/mgeradin/JuDyn.jl/blob/{commit}{path}#{line}",
    sitename="JuDyn.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
