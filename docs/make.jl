using Documenter
# using DocumenterLaTeX
# include("../src/JuDyn.jl")
# using .JUDYN
using JuDyn
push!(LOAD_PATH,"./src/")
push!(LOAD_PATH,"./src/MyAlgebra/src/")
push!(LOAD_PATH,"./src/set_elements/src/")
# push!(LOAD_PATH,"./src/solve/")
# makedocs(format = DocumenterLaTeX.LaTeX(),sitename="JuDyn")
makedocs(sitename="JuDyn")