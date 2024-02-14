using Documenter
# using DocumenterLaTeX
# include("../src/JuDyn.jl")
# using .JUDYN
using JuDyn
push!(LOAD_PATH,"./src/")
# makedocs(format = DocumenterLaTeX.LaTeX(),sitename="JuDyn")
# pages = ["Home" => "index.md", "Manual" => "intro.md"]
makedocs(
pages = [ "Home " => "index.md", "User Guide" => Any[ "Introduction" => "intro.md", 
"Modules"  => "modules.md"]], 
sitename="JuDyn")