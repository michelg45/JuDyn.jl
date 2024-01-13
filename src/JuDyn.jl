__precompile__()

module JuDyn

using Revise
using Reexport
using JLD
using JSON

include("MyAlgebra.jl")

@reexport using .MyAlgebra

include("FORDYN.jl")
@reexport using .FORDYN

include("Nodes.jl")
@reexport using .Nodes

include("Utils.jl")
@reexport using .Utils

include("LinearBeam.jl")
@reexport using .LinearBeam

include("SetElements.jl")
@reexport using .SetElements

include("BoundaryConditions.jl")
@reexport using .BoundaryConditions

include("Assemble.jl")
@reexport using .Assemble

include("InitialConditions.jl")
@reexport using .InitialConditions

include("Frames.jl")
@reexport using .Frames

include("BuildElements.jl")
@reexport using .BuildElements

include("InputFunctions.jl")
@reexport using .InputFunctions

include("Solve.jl")
@reexport using .Solve

end