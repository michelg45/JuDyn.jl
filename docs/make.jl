using Documenter
# using DocumenterLaTeX
# include("../src/JuDyn.jl")
# using .JUDYN
using JuDyn
push!(LOAD_PATH,"./src/")
# makedocs(format = DocumenterLaTeX.LaTeX(),sitename="JuDyn")
# pages = ["Home" => "index.md", "Manual" => "intro.md"]
makedocs(
pages = [ "Home " => "index.md", "User Guide" => Any[ "Introduction" => "intro.md", "Installation procedure" => "install.md", "Modules"  => "modules.md", "MyAlgebra" => "my_algebra.md", "Nodes" => "nodes.md","Set Elements" => "set_elements.md", "Element construction" => "build_elements.md"]], 
sitename="JuDyn")