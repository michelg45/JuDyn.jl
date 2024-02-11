using Documenter
# using DocumenterLaTeX
# include("../src/JuDyn.jl")
# using .JUDYN
using JuDyn
push!(LOAD_PATH,"./src/")
# makedocs(format = DocumenterLaTeX.LaTeX(),sitename="JuDyn")
# pages = ["Home" => "index.md", "Manual" => "intro.md"]
makedocs(
pages = Any["Home " => "index.md", "Introduction" => "intro.md"], 
sitename="JuDyn")